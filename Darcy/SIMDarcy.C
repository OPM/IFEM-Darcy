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
#include "CompatibleDarcy.h"

#include "AnaSol.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "Profiler.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMenums.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <cmath>
#include <strings.h>
#include "tinyxml2.h"


template<class Dim>
SIMDarcy<Dim>::SIMDarcy (Darcy& itg, unsigned char nf) :
  SIMMultiPatchModelGen<Dim>(nf), drc(itg), solVec(nullptr)
{
  this->initProblem();
}


template<class Dim>
SIMDarcy<Dim>::SIMDarcy (Darcy& itg, const CharVec& nf) :
  SIMMultiPatchModelGen<Dim>(nf), drc(itg), solVec(nullptr)
{
  this->initProblem();
}


template<class Dim>
void SIMDarcy<Dim>::initProblem ()
{
  drc.setOwnerSim(this);
  Dim::myProblem = &drc;
  Dim::myHeading = "Darcy solver";
  aCode[0] = aCode[1] = 0;

  if (drc.getOrder() > 0) // transient problem
    Dim::msgLevel = 1; // prints primary solution summary only
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
bool SIMDarcy<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"darcy"))
    return this->Dim::parse(elem);

  if (bool useCache = false; utl::getAttribute(elem,"cache",useCache)) {
    IFEM::cout << (useCache ? "\tEnabling" : "\tDisabling")
               <<" caching of element matrices."<< std::endl;
    drc.lCache(useCache);
  }

  bool gotMaterialData = false;
  DarcyMaterial defaultMaterial;
  if (mVec.empty()) mVec.reserve(1);

  for (const tinyxml2::XMLElement* child = elem->FirstChildElement();
       child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"materialdata")) {
      IFEM::cout <<"\tMaterial data with code "
                 << this->parseMaterialSet(child,mVec.size()) <<":\n";
      mVec.emplace_back(child);
    }
    else if (defaultMaterial.parse(child))
      gotMaterialData = true;
    else if (!strcasecmp(child->Value(),"anasol")) {
      IFEM::cout <<"  Parsing <anasol>"<< std::endl;
      Dim::mySol = new AnaSol(child, !this->mixedProblem());

      // Define the analytical boundary traction field
      if (int code = 0; utl::getAttribute(child,"code",code))
        if (code > 0 && Dim::mySol->getScalarSecSol())
        {
          this->setPropertyType(code,Property::NEUMANN);
          Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          aCode[1] = code;
        }
    }
    else if (!strcasecmp(child->Value(),"subiterations")) {
      IFEM::cout <<"  Parsing <subiterations>";
      utl::getAttribute(child,"tol",cycleTol);
      utl::getAttribute(child,"max",maxCycle);
      IFEM::cout <<"\n\ttol = "<< cycleTol <<" max = "<< maxCycle << std::endl;
    }
    else if (!Dim::myProblem->parse(child))
      this->Dim::parse(child);

  if (gotMaterialData && mVec.empty())
    mVec.push_back(std::move(defaultMaterial));

  bool ok = true;
  for (const DarcyMaterial& mat : mVec)
    if (mat.getViscosity() < 0.0)
      ok = false;
    else if (Matrix K; mat.getPermeability(&K))
      if (K.rows() != Dim::dimension || K.cols() != Dim::dimension)
      {
        std::cerr <<" *** SIMDarcy::parse(): Invalid permeability"
                  <<" matrix dimension ("<< K.rows() <<"x"<< K.cols() <<"),"
                  <<" should be ("<< Dim::dimension <<"x"<< Dim::dimension
                  <<")."<< std::endl;
        ok = false;
      }

  return ok;
}


template<class Dim>
bool SIMDarcy<Dim>::initNeumann (size_t propInd)
{
  if (const auto sit = Dim::myScalars.find(propInd);
      sit != Dim::myScalars.end())
    drc.setFlux(sit->second);

  else if (const auto tit = Dim::myTracs.find(propInd);
           tit != Dim::myTracs.end())
    drc.setFlux(tit->second);

  else if (const auto vit = Dim::myVectors.find(propInd);
           vit != Dim::myVectors.end())
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
  if (aCode[0] > 0)
    Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0)
    Dim::myVectors.erase(aCode[1]);
  aCode[0] = aCode[1] = 0;

  drc.setFlux(static_cast<RealFunc*>(nullptr));
  drc.setFlux(static_cast<VecFunc*>(nullptr));
  newTangent = true;
  mVec.clear();
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
                         Dim::opt.saveNorms ? &eNorm : nullptr);
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

  const int iType = tp.multiSteps() ? 0 : 1;
  const int iDump = iType == 0 ? tp.step/Dim::opt.saveInc : 1;

  if (newSolution && !solVec->empty())
  {
    // Write solution fields
    bool ok = this->writeGlvS(*solVec,iDump,nBlock,tp.time.t);

    // Write Darcy flux (velocity) vectors
    if (this->mixedProblem())
      ok &= this->writeGlvV(*solVec,"velocity",iDump,nBlock,110,Dim::nsd);
    else if (!Dim::opt.pSolOnly)
    {
      // Calculate Darcy flux vectors by projecting the secondary solution
      Matrix tmp;
      drc.set2ndFluxOnly(true);
      ok &= this->project(tmp,*solVec);
      drc.set2ndFluxOnly(false);
      ok &= this->writeGlvV(tmp,"velocity",iDump,nBlock,110,Dim::nsd);
    }
    if (!ok) return false;

    if (!Dim::opt.pSolOnly)
    {
      // Project the secondary solution onto the splines basis
      Vectors::iterator sit = proj.begin();
      for (const SIMoptions::ProjectionMap::value_type& pit : Dim::opt.project)
        if (solVec != &solution.front() &&
            !this->project(*sit,*solVec,pit.first))
          return false;
        else if (!this->writeGlvP(*sit,iDump,nBlock,100,pit.second.c_str()))
          return false;
        else
          ++sit;
    }

    // Write element norms
    if (Dim::opt.saveNorms)
      if (!this->writeGlvN(eNorm,iDump,nBlock))
        return false;
  }

  return this->writeGlvStep(iDump, iType == 1 ? iDump : tp.time.t, iType);
}


template<class Dim>
bool SIMDarcy<Dim>::init ()
{
  this->initSolution(this->getNoDOFs(), 1 + drc.getOrder());

  if (!solVec) solVec = &solution.front();
  this->registerField(this->mixedProblem() ? "solution" : "pressure", *solVec);

  this->initSystem(Dim::opt.solver);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);

  if (mVec.empty())
    mVec.emplace_back();
  if (mVec.size() == 1)
    drc.setMaterial(mVec.front());

  return true;
}


template<class Dim>
void SIMDarcy<Dim>::keepStep (const TimeStep& tp)
{
  if (Dim::msgLevel >= 0 && tp.multiSteps())
    this->printStep(tp.step, tp.time);

  newSolution = false;
}


template<class Dim>
bool SIMDarcy<Dim>::solveStep (const TimeStep& tp)
{
  if (Dim::msgLevel >= 0 && tp.multiSteps())
    this->printStep(tp.step, tp.time);

  if (tp.step == drc.getOrder())
    this->initLHSbuffers();

  const int printSol = maxCycle > -1 || tp.multiSteps() ? 0 : Dim::msgLevel;
  const int oldLevel = Dim::msgLevel;
  if (maxCycle > -1)
    Dim::msgLevel = 1; // suppress patch assembly loop log

  for (newSolution = false, tp.iter = 0; !newSolution; tp.iter++)
  {
    if (!this->setMode(SIM::DYNAMIC))
      return false;

    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);
    if (!this->assembleSystem(tp.time, solution, newTangent))
      return false;

    if (!this->solveSystem(solution.front()))
      return false;

    this->printSolutionSummary(solution.front(),printSol,nullptr,0);

    if (maxCycle < 0)
      newSolution = true;
    else
    {
      this->setMode(SIM::RHS_ONLY);
      this->updateDirichlet(tp.time.t,nullptr);
      this->assembleSystem(tp.time, solution);
      Vector linRes;
      this->extractLoadVec(linRes);
      double rConv = linRes.norm2();

      Dim::adm.cout <<"  cycle "<< tp.iter <<": Res = "<< rConv << std::endl;
      if (rConv < cycleTol)
        newSolution = true;
      else if (tp.iter >= maxCycle) {
        std::cerr <<" *** SIMDarcy::solveStep: Did not converge in "
                  << maxCycle <<" staggering cycles, bailing.."<< std::endl;
        break;
      }
    }

    newTangent = tp.step < drc.getOrder() || !drc.lCache();
  }

  if (maxCycle > -1)
  {
    this->printSolutionSummary(solution.front(),printSol,nullptr,0);
    Dim::msgLevel = oldLevel; // in case it is used by other SIM's
  }

  if (solVec == &solution.front()) // Not for adaptive simulations
  {
    // Project the secondary solution onto the splines basis
    Vectors::iterator sit = proj.begin();
    for (const SIMoptions::ProjectionMap::value_type& pit : Dim::opt.project)
      if (!this->project(*sit,*solVec,pit.first))
        return false;
  }

  return newSolution;
}


template<class Dim>
bool SIMDarcy<Dim>::advanceStep (TimeStep&)
{
  if (drc.advanceStep())
    this->pushSolution();

  return true;
}


template<class Dim>
void SIMDarcy<Dim>::printSolutionSummary (const Vector& solution,
                                          int printSol, const char*,
                                          std::streamsize outPrec)
{
  if (this->mixedProblem())
  {
    constexpr size_t nsd = Dim::dimension;

    // Compute and print solution norms for each component
    size_t iMax[Dim::dimension+1];
    double dMax[Dim::dimension+1];
    char D = 'P';
    double dNorm1 = this->solutionNorms(solution,dMax,iMax,1,'D');
    double dNorm2 = this->solutionNorms(solution,dMax+1,iMax+1,1,D);
    double dNormq = hypot(dNorm1,dNorm2);
    if (nsd == 3)
      dNormq = hypot(dNormq,this->solutionNorms(solution,dMax+2,iMax+2,1,++D));
    double dNormp = this->solutionNorms(solution,dMax+nsd,iMax+nsd,1,++D);

    int oldPrec = Dim::adm.cout.precision();
    if (outPrec > 0)
      Dim::adm.cout << std::setprecision(outPrec);

    Dim::adm.cout <<"\n>>> Primary solution summary <<<\n  L2-norm (q,p)  : "
                  << utl::trunc(dNormq) <<" "<< utl::trunc(dNormp);
    for (size_t d = 0; d < nsd; d++)
      if (utl::trunc(dMax[d]) != 0.0)
        Dim::adm.cout <<"\n  Max "<< char('X'+d) <<"-velocity : "
                      << dMax[d] <<" node "<< iMax[d];
    if (utl::trunc(dMax[nsd]) != 0.0)
      Dim::adm.cout <<"\n  Max pressure   : "
                    << dMax[nsd] <<" node "<< iMax[nsd];
    Dim::adm.cout << std::endl;

    if (outPrec > 0)
      Dim::adm.cout << std::setprecision(oldPrec);
  }
  else if (this->getNoFields() == 2)
  {
    // Compute and print solution norms
    size_t iMax[2];
    double dMax[2];
    double dNorm = this->solutionNorms(solution,dMax,iMax);

    int oldPrec = Dim::adm.cout.precision();
    if (outPrec > 0)
      Dim::adm.cout << std::setprecision(outPrec);

    Dim::adm.cout <<"  Primary solution summary: L2-norm      : "
                  << utl::trunc(dNorm)
                  <<"\n                            Max pressure : "
                  << utl::trunc(dMax[0]) <<" node "<< iMax[0]
                  <<"\n                       Max concentration : "
                  << utl::trunc(dMax[1]) <<" node "<< iMax[1] << std::endl;

    if (outPrec > 0)
      Dim::adm.cout << std::setprecision(oldPrec);
  }
  else
    this->SIMbase::printSolutionSummary(solution,printSol,"pressure",outPrec);
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


namespace
{
  //! \brief Prints some solution norms to log stream.
  void printSolNorms (const Vector& gNorm, bool compatible, size_t w)
  {
    if (compatible)
      if (double gn = gNorm[CompatibleDarcyNorm::H1_Qh]; utl::trunc(gn) != 0.0)
        IFEM::cout <<"\n  H1 norm |q^h| = a(q^h,q^h)^0.5"
                   << utl::adjustRight(w-32,"") << gn;
    if (double gn = gNorm[DarcyNorm::H1_Ph]; utl::trunc(gn) != 0.0)
      IFEM::cout << (compatible ?
                     "\n  L2 norm |p^h| = (p^h,p^h)^0.5 " :
                     "\n  H1 norm |p^h| = a(p^h,p^h)^0.5")
                 << utl::adjustRight(w-32,"") << gn;
    if (double gn = gNorm[DarcyNorm::EXT_ENERGY]; utl::trunc(gn) != 0.0)
      IFEM::cout <<"\n  External energy |(h,p^h)|^0.5"
                 << utl::adjustRight(w-31,"") << gn;
    if (double gn = gNorm[DarcyNorm::H1_Ch]; utl::trunc(gn) != 0.0)
      IFEM::cout <<"\n  H1 norm |c^h| = a(c^h,c^h)^0.5"
                 << utl::adjustRight(w-32,"") << gn;
  }

  //! \brief Prints some norms related to the exact solution to log stream.
  void printExactNorms (const Vector& gNorm, bool compatible, size_t w)
  {
    const double q   = compatible ? gNorm[CompatibleDarcyNorm::H1_Q] : 0.0;
    const double eqh = compatible ? gNorm[CompatibleDarcyNorm::H1_E_Qh] : 0.0;
    const double p   = gNorm[DarcyNorm::H1_P];
    const double eph = gNorm[DarcyNorm::H1_E_Ph];
    const double c   = gNorm[DarcyNorm::H1_C];
    const double ech = gNorm[DarcyNorm::H1_E_Ch];

    if (compatible && utl::trunc(q) != 0.0)
      IFEM::cout <<"\n  H1 norm |q| = a(q,q)^0.5"
                 << utl::adjustRight(w-26,"") << q;
    if (compatible && utl::trunc(eqh) != 0.0)
      IFEM::cout <<"\n  H1 norm |e| = a(e,e)^0.5, e=q-q^h"
                 << utl::adjustRight(w-35,"") << eqh;
    if (utl::trunc(p) != 0.0)
      IFEM::cout << (compatible ?
                     "\n  L2 norm |p| = (p,p)^0.5 " :
                     "\n  H1 norm |p| = a(p,p)^0.5")
                 << utl::adjustRight(w-26,"") << p;
    if (utl::trunc(eph) != 0.0)
      IFEM::cout << (compatible ?
                     "\n  L2 norm |e| = (e,e)^0.5, e=p-p^h " :
                     "\n  H1 norm |e| = a(e,e)^0.5, e=p-p^h")
                 << utl::adjustRight(w-35,"") << eph;
    if (utl::trunc(c) != 0.0)
      IFEM::cout <<"\n  H1 norm |c| = a(c,c)^0.5"
                 << utl::adjustRight(w-26,"") << c;
    if (utl::trunc(ech) != 0.0)
      IFEM::cout <<"\n  H1 norm |e| = a(e,e)^0.5, e=c-c^h"
                 << utl::adjustRight(w-35,"") << ech;
    if (compatible && fabs(q) > 1.0e-16)
      IFEM::cout <<"\n  Exact relative error (in % of q)"
                 << utl::adjustRight(w-34,"") << 100.0*eqh / q;
    if (double pc = hypot(p,c); pc > 1.0e-16)
      IFEM::cout <<"\n  Exact relative error ("<< (compatible ? "in % of p":"%")
                 <<")"<< utl::adjustRight(compatible ? w-34 : w-26, "")
                 << 100.0*hypot(eph,ech) / pc;
  }

  //! \brief Prints a norm group to the log stream.
  void printNormGroup (const Vector& rNorm, const Vector& fNorm,
                       const std::string& name, bool noAnaSol, size_t w)
  {
    IFEM::cout <<"\nError estimates based on >>> "<< name <<" <<<";
    if (double rn = rNorm[DarcyNorm::H1_Pr_Ph]; utl::trunc(rn) != 0.0)
      IFEM::cout <<"\n  H1 norm |p^r-p^h|"<< utl::adjustRight(w-19,"") << rn;
    if (double rn = rNorm[DarcyNorm::H1_Cr_Ch]; utl::trunc(rn) != 0.0)
      IFEM::cout <<"\n  H1 norm |c^r-c^h|"<< utl::adjustRight(w-19,"") << rn;

    if (noAnaSol)
      return;

    if (double rn = rNorm[DarcyNorm::H1_E_Pr]; utl::trunc(rn) != 0.0)
      IFEM::cout <<"\n  H1 norm |p^r-p|"<< utl::adjustRight(w-17,"") << rn;
    if (double rn = rNorm[DarcyNorm::H1_E_Cr]; utl::trunc(rn) != 0.0)
      IFEM::cout <<"\n  H1 norm |c^r-c|"<< utl::adjustRight(w-17,"") << rn;

    if (double fn = fNorm[DarcyNorm::H1_E_Ph]; fabs(fn) > 1.0e-16)
      IFEM::cout <<"\n  Effectivity index eta^p"<< utl::adjustRight(w-25,"")
                 << rNorm[DarcyNorm::H1_Pr_Ph] / fn;
    if (double fn = fNorm[DarcyNorm::H1_E_Ch]; fabs(fn) > 1.0e-16)
      IFEM::cout <<"\n  Effectivity index eta^c"<< utl::adjustRight(w-25,"")
                 << rNorm[DarcyNorm::H1_Cr_Ch] / fn;
    if (double fn = fNorm[DarcyNorm::TOTAL_NORM_E]; fabs(fn) > 1.0e-16)
      IFEM::cout <<"\n  Effectivity index eta^tot"<< utl::adjustRight(w-27,"")
                 << rNorm[DarcyNorm::TOTAL_NORM_REC] / fn;
  }
}


template<class Dim>
void SIMDarcy<Dim>::printFinalNorms (const TimeStep& tp)
{
  // Don't print final norms with adaptive simulations
  if (solVec != &solution.front())
    return;

  if (!this->setMode(SIM::RECOVERY))
    return;

  // Evaluate solution norms and print to terminal
  Vectors gNorm;
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (this->solutionNorms(tp.time,{solution.front()},proj,gNorm))
    this->printNorms(gNorm,36);
}


template<class Dim>
void SIMDarcy<Dim>::printNorms (const Vectors& gNorm, size_t w) const
{
  if (gNorm.empty())
    return;

  IFEM::cout <<"\n>>> Norm summary <<<";
  printSolNorms(gNorm.front(),this->mixedProblem(),w);

  if (Dim::mySol)
    printExactNorms(gNorm.front(),this->mixedProblem(),w);

  Vectors::const_iterator git = gNorm.begin();
  for (const SIMoptions::ProjectionMap::value_type& pit : Dim::opt.project)
    if (++git != gNorm.end())
      printNormGroup(*git, gNorm.front(), pit.second, !Dim::mySol, w);

  IFEM::cout << std::endl;
}


/*!
  This method is overridden to resolve inhomogeneous boundary condition fields,
  in case they are derived from the analytical solution.
*/

template<class Dim>
bool SIMDarcy<Dim>::preprocessA ()
{
  proj.resize(Dim::opt.project.size());
  if (!Dim::mySol) return true;

  // Define analytical boundary condition fields
  for (Property& prop : Dim::myProps)
    if (prop.pcode == Property::DIRICHLET_ANASOL)
    {
      if (!Dim::mySol->getScalarSol())
        prop.pcode = Property::UNDEFINED;
      else if (aCode[0] == abs(prop.pindx))
        prop.pcode = Property::DIRICHLET_INHOM;
      else if (aCode[0] == 0)
      {
        aCode[0] = abs(prop.pindx);
        Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
        prop.pcode = Property::DIRICHLET_INHOM;
      }
      else
        prop.pcode = Property::UNDEFINED;
    }
    else if (prop.pcode == Property::NEUMANN_ANASOL)
    {
      if (!Dim::mySol->getScalarSecSol())
        prop.pcode = Property::UNDEFINED;
      else if (aCode[1] == prop.pindx)
        prop.pcode = Property::NEUMANN;
      else if (aCode[1] == 0)
      {
        aCode[1] = prop.pindx;
        Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
        prop.pcode = Property::NEUMANN;
      }
      else
        prop.pcode = Property::UNDEFINED;
    }

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::preprocessB ()
{
  if (this->getNoConstraints() == 0 && !drc.extEner)
    drc.extEner = 'y';
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


template class SIMDarcy<SIM1D>;
template class SIMDarcy<SIM2D>;
template class SIMDarcy<SIM3D>;
