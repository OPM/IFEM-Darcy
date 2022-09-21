// $Id$
//==============================================================================
//!
//! \file SIMDarcy.C
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Simulation driver for Isogeometric FE analysis of Darcy Flow.
//!
//==============================================================================

#include "SIMDarcy.h"

#include "Darcy.h"
#include "DarcySolutions.h"

#include "AnaSol.h"
#include "DataExporter.h"
#include "ExprFunctions.h"
#include "Functions.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Profiler.h"
#include "Property.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMenums.h"
#include "TimeDomain.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <cmath>
#include <cstdlib>
#include <strings.h>
#include "tinyxml.h"


template<class Dim>
SIMDarcy<Dim>::SIMDarcy (Darcy& itg, unsigned char nf) :
  SIMDarcy<Dim>(itg, std::vector<unsigned char>{nf})
{
}


template<class Dim>
SIMDarcy<Dim>::SIMDarcy (Darcy& itg, const std::vector<unsigned char>& nf) :
  SIMMultiPatchModelGen<Dim>(nf), drc(itg), solVec(nullptr)
{
  Dim::myProblem = &drc;
  Dim::myHeading = "Darcy solver";
  aCode[0] = aCode[1] = 0;
}


template<class Dim>
SIMDarcy<Dim>::~SIMDarcy ()
{
  Dim::myProblem = nullptr;
  Dim::myInts.clear();
  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
}


template<class Dim>
bool SIMDarcy<Dim>::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"darcy"))
    return this->Dim::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    const char* value = nullptr;
    if ((value = utl::getValue(child,"bodyforce")))
      drc.setBodyForce(new VecFuncExpr(value));
    else if ((value = utl::getValue(child,"gravity")))
      drc.setGravity(atof(value));
    else if (!strcasecmp(child->Value(), "materialdata")) {
      int code = this->parseMaterialSet(child,mVec.size());
      mVec.resize(mVec.size()+1);
      IFEM::cout << "\tMaterial data with code " << code <<":\n";
      if (!mVec.back().parse(child))
        mVec.pop_back();
    } else if (DarcyMaterial::handlesTag(child->Value())) {
      if (mVec.empty())
        mVec.resize(1);
      mVec.back().parse(child);
    } else if (!strcasecmp(child->Value(),"source") || !strcasecmp(child->Value(),"source_c")) {
      bool isC = strcasecmp(child->Value(),"source") != 0;
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\tSource function" << (isC ? " (concentration):" : ":");
      RealFunc* src = nullptr;
      const char* input = isC ? utl::getValue(child, "source_c")
                              : utl::getValue(child, "source");
      if (type == "expression") {
        if (input)
          IFEM::cout << " " << input << std::endl;
        src = new EvalFunction(input);
      }
      else if (type == "diracsum") {
        double tol = 1e-2;
        utl::getAttribute(child, "pointTol", tol);
        if (input) {
          IFEM::cout << " DiracSum";
          DiracSum* f = new DiracSum(tol, Dim::dimension);
          if (f->parse(input))
            src = f;
          else
            delete f;
        }
      }
      else if (type == "elementsum") {
        if (input) {
          if (!this->createFEMmodel())
            continue;
          IFEM::cout << " ElementSum";
          ElementSum* f = new ElementSum(Dim::dimension);

          if (f->parse(input, *this))
            src = f;
          else
            delete f;
        }
      } else
        IFEM::cout <<"(none)"<< std::endl;

      if (src) {
        if (isC)
          drc.setCSource(src);
        else
          drc.setSource(src);
      }
    }
    else if (!strcasecmp(child->Value(),"anasol")) {
      std::string type;
      utl::getAttribute(child,"type",type);
      if (type == "Lshape") {
        Dim::mySol = new AnaSol(new LshapeDarcy(), new LshapeDarcyVelocity());
        IFEM::cout <<"\tAnalytical solution: Lshape"<< std::endl;
      }
      else {
        Dim::mySol = new AnaSol(child);
        IFEM::cout <<"\tAnalytical solution: expression"<< std::endl;
      }
      // Define the analytical boundary traction field
      int code = 0;
      if (utl::getAttribute(child,"code",code) && code > 0) {
        if (code > 0 && Dim::mySol->getScalarSecSol())
        {
          this->setPropertyType(code,Property::NEUMANN);
          Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          aCode[1] = code;
        }
      }
    }
    else if (!strcasecmp(child->Value(),"reactions"))
      drc.extEner = 'R';
    else if (!strcasecmp(child->Value(),"subiterations")) {
      IFEM::cout << "\tUsing sub-iterations";
      utl::getAttribute(child,"tol",cycleTol);
      utl::getAttribute(child,"max",maxCycle);
      IFEM::cout <<"\n\t\ttol = "<< cycleTol;
      IFEM::cout <<"\n\t\tmax = "<< maxCycle;
      IFEM::cout << std::endl;
    }
    else
      this->Dim::parse(child);
  }

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::initNeumann (size_t propInd)
{
  const auto sit = Dim::myScalars.find(propInd);
  const auto tit = Dim::myTracs.find(propInd);
  const auto vit = Dim::myVectors.find(propInd);

  if (sit != Dim::myScalars.end())
    drc.setFlux(sit->second);
  else if (tit != Dim::myTracs.end())
    drc.setFlux(tit->second);
  else if (vit != Dim::myVectors.end())
    drc.setFlux(vit->second);
  else
    return false;

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size())
    return false;

  drc.setMaterial(mVec[propInd]);

  return true;
}


template<class Dim>
void SIMDarcy<Dim>::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  aCode[0] = aCode[1] = 0;

  drc.setFlux((RealFunc*)nullptr);
  drc.setFlux((VecFunc*)nullptr);
  this->Dim::clearProperties();
}


template<class Dim>
void SIMDarcy<Dim>::registerFields (DataExporter& exporter)
{
  int results = DataExporter::PRIMARY;

  if (!Dim::opt.pSolOnly)
    results |= DataExporter::SECONDARY;

  if (Dim::opt.saveNorms)
    results |= DataExporter::NORMS;

  exporter.registerField("u", "primary", DataExporter::SIM, results);
  exporter.setFieldValue("u", this, solVec,
                       Dim::opt.project.empty() ? nullptr : &proj,
                       (results & DataExporter::NORMS) ? &eNorm : nullptr);
}


template<class Dim>
bool SIMDarcy<Dim>::saveModel (char* fileName, int& geoBlk, int& nBlock)
{
  if (Dim::opt.format < 0)
    return true;

  nBlock = 0;
  return this->writeGlvG(geoBlk,fileName);
}


template<class Dim>
bool SIMDarcy<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  if (Dim::opt.format < 0 || (tp.step % Dim::opt.saveInc) > 0)
    return true;

  int iDump = tp.step/Dim::opt.saveInc + (drc.getOrder() == 0 ? 1 : 0);

  // Write solution fields
  if (!this->writeGlvS(*solVec,iDump,nBlock,tp.time.t))
    return false;

  if (!this->doProjection())
    return false;

  if (!solVec->empty())  {
    if (!Dim::opt.pSolOnly) {
      Matrix tmp;
      if (!this->project(tmp,*solVec))
        return false;

      if (!this->writeGlvV(tmp,"velocity",iDump,nBlock,110,Dim::nsd))
        return false;

      size_t i = 0;
      for (auto& pit : Dim::opt.project)
        if (!this->writeGlvP(proj[i++],iDump,nBlock,100,pit.second.c_str()))
          return false;
    }
  }

  // Write element norms
  if (Dim::opt.saveNorms)
    if (!this->writeGlvN(eNorm,iDump,nBlock))
      return false;

  double param2 = drc.getOrder() == 0 ? iDump : tp.time.t;
  return this->writeGlvStep(iDump,param2,drc.getOrder() == 0 ? 1 : 0);
}


template<class Dim>
bool SIMDarcy<Dim>::init ()
{
  this->initSolution(this->getNoDOFs(), 1 + drc.getOrder());
  if (!this->solVec)
    this->solVec = &this->solution.front();

  this->initSystem(Dim::opt.solver);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);

  if (mVec.empty())
    mVec.push_back(DarcyMaterial());
  if (mVec.size() == 1)
    drc.setMaterial(mVec.front());

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::solveStep (const TimeStep& tp)
{
  if (Dim::msgLevel >= 0 && tp.multiSteps())
    IFEM::cout <<"\n  step = "<< tp.step
               <<"  time = "<< tp.time.t << std::endl;

  bool conv = false;
  tp.iter = 0;
  while (!conv)
  {
    if (!this->setMode(SIM::DYNAMIC))
      return false;

    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);
    if (!this->assembleSystem(tp.time, solution))
      return false;

    if (!this->solveSystem(solution.front(),maxCycle > -1 ? 0 : Dim::msgLevel-1,"pressure    "))
      return false;

    if (maxCycle > -1) {
      this->setMode(SIM::RHS_ONLY);
      this->updateDirichlet(tp.time.t,nullptr);
      this->assembleSystem(tp.time, solution);
      Vector linRes;
      this->extractLoadVec(linRes);
      double rConv = linRes.norm2();

      IFEM::cout <<"  cycle "<< tp.iter <<": Res = "<< rConv << std::endl;
      if (rConv < cycleTol)
        conv = true;

      if (tp.iter >= maxCycle) {
        std::cerr <<" *** SIMDarcy::solveStep: Did not converge in "
                 << maxCycle <<" staggering cycles, bailing.."<< std::endl;
        return false;
      }
    } else
      conv = true;

    ++tp.iter;
  }

  if (maxCycle > -1)
    this->printSolutionSummary(solution.front(), 0, nullptr, 0);

  if (!tp.multiSteps() && adNorm == DCY::NO_ADAP && !this->doProjection())
    return false;

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::advanceStep (TimeStep&)
{
  if (drc.getOrder() > 0) {
    drc.advanceStep();
    this->pushSolution();
  }
  return true;
}


template<class Dim>
void SIMDarcy<Dim>::printSolutionSummary (const Vector& solution,
                                          int printSol, const char*,
                                          std::streamsize outPrec)
{
  const size_t nf = this->getNoFields();
  if (nf > 1) {
    // Compute and print solution norms
    size_t iMax[nf];
    double dMax[nf];
    double dNorm;
    drc.getSolutionNorms(*this,solution,dNorm,dMax,iMax);

    int oldPrec = this->adm.cout.precision();
    if (outPrec > 0)
      this->adm.cout << std::setprecision(outPrec);

    this->adm.cout <<"  Primary solution summary: L2-norm         : ";
    this->adm.cout << utl::trunc(dNorm);
    this->adm.cout <<"\n                               Max pressure : ";
    this->adm.cout << dMax[0] <<" node "<< iMax[0];
    this->adm.cout <<"\n                          Max concentration : ";
    this->adm.cout << dMax[1] <<" node "<< iMax[1] << "\n";

    this->adm.cout << std::setprecision(oldPrec);
  } else {
    this->SIMbase::printSolutionSummary(solution, printSol,
                                        "pressure    ", outPrec);
  }
}


template<class Dim>
bool SIMDarcy<Dim>::solveSystem (Vector& solution, int printSol,
                                 double* rCond, const char* compName,
                                 size_t idxRHS)
{
  if (!this->Dim::solveSystem(solution,printSol,rCond,compName,idxRHS))
    return false;
  else if (idxRHS > 0 || !this->haveReactions() || drc.extEner != 'R')
    return true;

  // Assemble the reaction forces. Strictly, we only need to assemble those
  // elements that have nodes on the Dirichlet boundaries, but...
  return this->assembleForces({solution},0.0,&myReact);
}


template<class Dim>
SIM::ConvStatus SIMDarcy<Dim>::solveIteration (TimeStep& tp)
{
  return this->solveStep(tp) ? SIM::CONVERGED : SIM::FAILURE;
}


template<class Dim>
void SIMDarcy<Dim>::printSolNorms (const Vector& gNorm,
                                   size_t w) const
{
  IFEM::cout << "\n  H1 norm |p^h| = a(p^h,p^h)^0.5"
             << utl::adjustRight(w-32,"") << gNorm[DarcyNorm::H1_Ph];
  if (utl::trunc(gNorm[DarcyNorm::EXT_ENERGY]) != 0.0)
    IFEM::cout << "\n  External energy |(h,p^h)|^0.5"
                << utl::adjustRight(w-31,"") << gNorm[DarcyNorm::EXT_ENERGY];
  if (this->getNoFields() > 1)
    IFEM::cout << "\n  H1 norm |c^h| = a(c^h,c^h)^0.5"
               << utl::adjustRight(w-32,"") << gNorm[DarcyNorm::H1_Ch];
}


template<class Dim>
void SIMDarcy<Dim>::printFinalNorms (const TimeStep& tp)
{
  // Don't print final norms with adaptive simulations
  if (solVec != &solution.front())
    return;

  if (!this->setMode(SIM::RECOVERY))
    return;

  if (!this->doProjection()) {
    std::cerr << "*** Error during solution projection." << std::endl;
    return;
  }

  if (!Dim::opt.project.empty())
  {
    // Project the secondary solution onto the splines basis
    size_t j = 0;
    for (auto& pit : Dim::opt.project)
      if (!this->project(proj[j++],solution.front(),pit.first))
        return;
  }

  // Evaluate solution norms
  Vectors gNorm;
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (!this->solutionNorms(tp.time,Vectors(1, solution.front()),proj,gNorm))
    return;

  // Print global norm summary to console
  this->printNorms(gNorm, 36);

}

template<class Dim>
void SIMDarcy<Dim>::printNorms (const Vectors& gNorm, size_t w) const
{
  if (gNorm.empty()) return;

  IFEM::cout << "\n>>> Norm summary <<<";
  this->printSolNorms(gNorm.front(),w);

  if (Dim::mySol)
    this->printExactNorms(gNorm.front(),w);

  size_t j = 0;
  for (const auto& prj : this->opt.project)
    if (++j < gNorm.size())
      this->printNormGroup(gNorm[j],gNorm[0],prj.second);

  IFEM::cout << std::endl;
}

template<class Dim>
void SIMDarcy<Dim>::printExactNorms (const Vector& gNorm,
                                     size_t w) const
{
  if (!Dim::mySol)
    return;

  IFEM::cout << "\n  H1 norm |p| = a(p,p)^0.5"
             << utl::adjustRight(w-26,"") << gNorm[DarcyNorm::H1_P];
  IFEM::cout << "\n  H1 norm |e| = a(e,e)^0.5, e=p-p^h"
             << utl::adjustRight(w-35,"") << gNorm[DarcyNorm::H1_E_Ph];
  if (this->getNoFields() > 1) {
    IFEM::cout << "\n  H1 norm |c| = a(c,c)^0.5"
               << utl::adjustRight(w-26,"") << gNorm[DarcyNorm::H1_C];
    IFEM::cout << "\n  H1 norm |e| = a(e,e)^0.5, e=c-c^h"
               << utl::adjustRight(w-35,"") << gNorm[DarcyNorm::H1_E_Ch];
    IFEM::cout << "\n  Exact relative error (%)"
               << utl::adjustRight(w-26,"") << hypot(gNorm[DarcyNorm::H1_E_Ph],gNorm[DarcyNorm::H1_E_Ch])*100.0 /
                                               hypot(gNorm[DarcyNorm::H1_P],gNorm[DarcyNorm::H1_C]);
  } else
    IFEM::cout << "\n  Exact relative error (%)"
               << utl::adjustRight(w-26,"") << gNorm[DarcyNorm::H1_E_Ph]*100.0 /
                                               gNorm[DarcyNorm::H1_P];
}


template<class Dim>
void SIMDarcy<Dim>::printNormGroup (const Vector& rNorm,
                                    const Vector& fNorm,
                                    const std::string& name) const
{
  IFEM::cout << "\nError estimates based on >>> " << name << " <<<";
  size_t w = 36;
  if (name == "Pure residuals")
    ; // TODO
  else {
    IFEM::cout << "\n  H1 norm |p^r-p^h|"
               << utl::adjustRight(w-19,"") << rNorm[DarcyNorm::H1_Pr_Ph];
    if (this->getNoFields() > 1)
      IFEM::cout << "\n  H1 norm |c^r-c^h|"
                 << utl::adjustRight(w-19,"") << rNorm[DarcyNorm::H1_Cr_Ch];
    if (Dim::mySol) {
      IFEM::cout << "\n  H1 norm |p^r-p|"
                 << utl::adjustRight(w-17,"") << rNorm[DarcyNorm::H1_E_Pr];
      if (this->getNoFields() > 1)
        IFEM::cout << "\n  H1 norm |c^r-c|"
                   << utl::adjustRight(w-17,"") << rNorm[DarcyNorm::H1_E_Cr];

      IFEM::cout << "\n  Effectivity index eta^p"
                 << utl::adjustRight(w-25,"")
                 << rNorm[DarcyNorm::H1_Pr_Ph] / fNorm[DarcyNorm::H1_E_Ph];
      if (this->getNoFields() > 1)
        IFEM::cout << "\n  Effectivity index eta^c"
                   << utl::adjustRight(w-25,"")
                   << rNorm[DarcyNorm::H1_Cr_Ch] / fNorm[DarcyNorm::H1_E_Ch]
                   << "\n  Effectivity index eta^tot"
                   << utl::adjustRight(w-27,"")
                   << rNorm[DarcyNorm::TOTAL_NORM_REC] / fNorm[DarcyNorm::TOTAL_NORM_E];
    }
  }
}


template<class Dim>
void SIMDarcy<Dim>::preprocessA ()
{
  proj.resize(Dim::opt.project.size());
  if (!Dim::mySol) return;

  // Define analytical boundary condition fields
  PropertyVec::iterator p;
  for (p = Dim::myProps.begin(); p != Dim::myProps.end(); ++p)
    if (p->pcode == Property::DIRICHLET_ANASOL)
    {
      if (!Dim::mySol->getScalarSol())
        p->pcode = Property::UNDEFINED;
      else if (aCode[0] == abs(p->pindx))
        p->pcode = Property::DIRICHLET_INHOM;
      else if (aCode[0] == 0)
      {
        aCode[0] = abs(p->pindx);
        Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
        p->pcode = Property::DIRICHLET_INHOM;
      }
      else
        p->pcode = Property::UNDEFINED;
    }
    else if (p->pcode == Property::NEUMANN_ANASOL)
    {
      if (!Dim::mySol->getScalarSecSol())
        p->pcode = Property::UNDEFINED;
      else if (aCode[1] == p->pindx)
        p->pcode = Property::NEUMANN;
      else if (aCode[1] == 0)
      {
        aCode[1] = p->pindx;
        Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
        p->pcode = Property::NEUMANN;
      }
      else
        p->pcode = Property::UNDEFINED;
    }
}


template<class Dim>
bool SIMDarcy<Dim>::preprocessB ()
{
  if (this->getNoConstraints() == 0 && !drc.extEner)
    drc.extEner = 'y';
  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::doProjection()
{
  if (!Dim::opt.project.empty() && proj.front().empty())
  {
    // Project the secondary solution onto the splines basis
    size_t j = 0;
    for (auto& pit : Dim::opt.project)
      if (!this->project(proj[j++],solution.front(),pit.first))
        return false;
  }

  return true;
}


template<class Dim>
double SIMDarcy<Dim>::getReferenceNorm (const Vectors& gNorm,
                                        size_t adaptor) const
{
  if (adaptor == 0) {
    if (adNorm == DCY::PRESSURE_H1)
      return gNorm[0][DarcyNorm::H1_P];
    else if (adNorm == DCY::CONCENTRATION_H1)
      return gNorm[0][DarcyNorm::H1_C];
    else
      return hypot(gNorm[0][DarcyNorm::H1_P], gNorm[0][DarcyNorm::H1_C]);
  }

  if (this->haveAnaSol()) {
    if (adNorm == DCY::RECOVERY_PRESSURE)
      return gNorm[0][DarcyNorm::H1_P];
    else if (adNorm == DCY::RECOVERY_CONCENTRATION)
      return gNorm[0][DarcyNorm::H1_C];
    else
      return hypot(gNorm[0][DarcyNorm::H1_P], gNorm[0][DarcyNorm::H1_C]);
  }

  return hypot(hypot(gNorm[0][DarcyNorm::H1_Ph], gNorm[0][DarcyNorm::H1_Ch]),
               hypot(gNorm[1][DarcyNorm::H1_Pr_Ph], gNorm[1][DarcyNorm::H1_Cr_Ch]));
}


template<class Dim>
double SIMDarcy<Dim>::getEffectivityIndex (const Vectors& gNorm,
                                           size_t adaptor,
                                           size_t inorm) const
{
  if (adNorm == DCY::RECOVERY_PRESSURE)
    return gNorm[adaptor](inorm) / gNorm[0][DarcyNorm::H1_E_Ph];

  if (adNorm == DCY::RECOVERY_CONCENTRATION)
    return gNorm[adaptor](inorm) / gNorm[0][DarcyNorm::H1_E_Ch];

  return gNorm[adaptor](inorm) / gNorm[0][DarcyNorm::TOTAL_NORM_E];
}


template<class Dim>
int SolverConfigurator<SIMDarcy<Dim>>::
setup (SIMDarcy<Dim>& darcy,
       const typename SIMDarcy<Dim>::SetupProps& props,
       char* infile)
{
  utl::profiler->start("Model input");

  if (!darcy.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (!darcy.preprocess())
    return 2;

  darcy.init();

  return 0;
}


#define INSTANCE(T) \
  template class SIMDarcy<T>; \
  template struct SolverConfigurator<SIMDarcy<T>>;

INSTANCE(SIM1D)
INSTANCE(SIM2D)
INSTANCE(SIM3D)
