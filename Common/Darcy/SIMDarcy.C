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

#include "DarcySolutions.h"

#include "AlgEqSystem.h"
#include "AnaSol.h"
#include "DataExporter.h"
#include "ExprFunctions.h"
#include "Functions.h"
#include "IFEM.h"
#include "Profiler.h"
#include "ReactionsOnly.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"

#include "tinyxml.h"


template<class Dim>
SIMDarcy<Dim>::SIMDarcy (int torder) :
  Dim(1), drc(Dim::dimension, torder), solVec(nullptr)
{
  Dim::myProblem = &drc;
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
    if ((value = utl::getValue(child,"permvalues"))) {
      IFEM::cout <<"\tPermeability: " << value << std::endl;
      drc.setPermValues(new VecFuncExpr(value));
    } else if ((value = utl::getValue(child,"permeability"))) {
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\tPermeability";
      drc.setPermField(utl::parseRealFunc(value,type));
      IFEM::cout << std::endl;
    }
    else if ((value = utl::getValue(child,"bodyforce")))
      drc.setBodyForce(new VecFuncExpr(value));
    else if ((value = utl::getValue(child,"gravity")))
      drc.setGravity(atof(value));
    else if (!strcasecmp(child->Value(),"source")) {
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\tSource function:";
      if (type == "expression" && child->FirstChild()) {
        IFEM::cout << " " << child->FirstChild()->Value() << std::endl;
        drc.setSource(new EvalFunction(child->FirstChild()->Value()));
      }
      else if (type == "diracsum") {
        double tol = 1e-2;
        const char* input = utl::getValue(child, "source");
        utl::getAttribute(child, "pointTol", tol);
        if (input) {
          IFEM::cout << " DiracSum";
          DiracSum* f = new DiracSum(tol, Dim::dimension);
          if (f->parse(input))
            drc.setSource(f);
          else
            delete f;
        }
      }
      else
        IFEM::cout <<"(none)"<< std::endl;
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
        std::cout <<"\tAnalytical solution: expression"<< std::endl;
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
    else
      this->Dim::parse(child);
  }

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::initNeumann (size_t propInd)
{
  typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
  typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

  if (sit != Dim::myScalars.end())
    drc.setFlux(sit->second);
  else if (vit != Dim::myVectors.end())
    drc.setFlux(vit->second);
  else
    return false;

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
  if (Dim::opt.format < 0)
    return true;

  int iDump = tp.step/Dim::opt.saveInc + (drc.getOrder() == 0 ? 1 : 0);

  // Write solution fields
  if (!this->writeGlvS(*solVec,iDump,nBlock,tp.time.t))
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
void SIMDarcy<Dim>::init ()
{
this->initSolution(this->getNoDOFs(), 1 + drc.getOrder());
this->solVec = &this->solution.front();

this->initSystem(Dim::opt.solver);
this->setQuadratureRule(Dim::opt.nGauss[0],true);
}


template<class Dim>
bool SIMDarcy<Dim>::solveStep (TimeStep& tp)
{
  if (Dim::msgLevel >= 0 && tp.multiSteps())
    IFEM::cout <<"\n  step = "<< tp.step
               <<"  time = "<< tp.time.t << std::endl;

  if (!this->setMode(SIM::DYNAMIC))
    return false;

  if (!this->assembleSystem(tp.time, solution))
    return false;

  if (!this->solveSystem(solution.front(),Dim::msgLevel-1,"pressure    "))
    return false;

  if (!this->setMode(SIM::RECOVERY))
    return false;

  if (!Dim::opt.project.empty() && !tp.multiSteps())
  {
    // Project the secondary solution onto the splines basis
    size_t j = 0;
    for (auto& pit : Dim::opt.project)
      if (!this->project(proj[j++],solution.front(),pit.first))
        return false;

    IFEM::cout << std::endl;
  }

  return true;
}


template<class Dim>
bool SIMDarcy<Dim>::advanceStep (TimeStep&)
{
  if (drc.getOrder() > 0)
    this->pushSolution();
  return true;
}


template<class Dim>
void SIMDarcy<Dim>::printSolutionSummary (const Vector& solution,
                                          int printSol, const char*,
                                          std::streamsize outPrec)
{
  this->SIMbase::printSolutionSummary(solution, printSol,
                                      "pressure    ", outPrec);
}


template<class Dim>
bool SIMDarcy<Dim>::solveSystem (Vector& solution, int printSol,
                                 double* rCond, const char* compName,
                                 bool newLHS, size_t idxRHS)
{
  if (!this->Dim::solveSystem(solution,printSol,rCond,compName,newLHS,idxRHS))
    return false;
  else if (idxRHS > 0 || !this->haveReactions() || drc.extEner != 'R')
    return true;

  // Assemble the reaction forces. Strictly, we only need to assemble those
  // elements that have nodes on the Dirichlet boundaries, but...
  drc.setReactionIntegral(new ReactionsOnly(myReact,Dim::mySam,Dim::adm));
  AlgEqSystem* tmpEqSys = Dim::myEqSys;
  Dim::myEqSys = nullptr;
  bool ok = this->setMode(SIM::RHS_ONLY) && this->assembleSystem({solution});
  Dim::myEqSys = tmpEqSys;
  drc.setReactionIntegral(nullptr);

  return ok;
}


template<class Dim>
void SIMDarcy<Dim>::printFinalNorms (const TimeStep& tp)
{
  // Evaluate solution norms
  Vectors gNorm;
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (tp.multiSteps()) {
    if (!this->solutionNorms(tp.time,solution,proj,gNorm))
      return;
  } else {
    if (!this->solutionNorms(*solVec,proj,gNorm))
      return;
  }

  // Print global norm summary to console
  this->printNorms(gNorm);
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
