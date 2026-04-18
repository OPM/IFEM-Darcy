// $Id$
//==============================================================================
//!
//! \file SIMDarcyAdvection.C
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of Darcy advection.
//!
//==============================================================================

#include "SIMDarcyAdvection.h"

#include "DarcyAdvection.h"
#include "DarcySolutions.h"

#include "DataExporter.h"
#include "IFEM.h"
#include "Profiler.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMenums.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <strings.h>
#include "tinyxml2.h"


template<class Dim>
SIMDarcyAdvection<Dim>::SIMDarcyAdvection (DarcyAdvection& itg) :
  Dim(1), drc(itg)
{
  Dim::myProblem = &drc;
  Dim::myHeading = "Darcy advection solver";
  Dim::msgLevel  = 1;
}


template<class Dim>
SIMDarcyAdvection<Dim>::~SIMDarcyAdvection ()
{
  Dim::myProblem = nullptr;
  Dim::myInts.clear();
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"darcyadvection"))
    return this->Dim::parse(elem);

  if (bool useCache = false; utl::getAttribute(elem,"cache",useCache)) {
    IFEM::cout << (useCache ? "\tEnabling" : "\tDisabling")
               <<" caching of element matrices."<< std::endl;
    drc.lCache(useCache);
  }

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(), "materialdata")) {
      int code = this->parseMaterialSet(child,mVec.size());
      mVec.resize(mVec.size()+1);
      IFEM::cout << "\tMaterial data with code " << code <<":\n";
      if (!mVec.back().parse(child))
        mVec.pop_back();
    }
    else if (DarcyMaterial::handlesTag(child->Value())) {
      if (mVec.empty())
        mVec.resize(1);
      mVec.back().parse(child);
    }
    else if (!Dim::myProblem->parse(child))
      this->Dim::parse(child);

  return true;
}


template<class Dim>
void SIMDarcyAdvection<Dim>::registerFields (DataExporter& exporter)
{
  int results = DataExporter::PRIMARY;

  if (!Dim::opt.pSolOnly)
    results |= DataExporter::SECONDARY;

  if (Dim::opt.saveNorms)
    results |= DataExporter::NORMS;

  exporter.registerField("c", "primary", DataExporter::SIM, results);
  exporter.setFieldValue("c", this, &solution.front());
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  if (Dim::opt.format < 0 || (tp.step % Dim::opt.saveInc) > 0)
    return true;

  int iDump = tp.step/Dim::opt.saveInc;

  // Write solution fields, step info is written by SIMDarcy
  return this->writeGlvS(solution.front(),iDump,nBlock,tp.time.t,"c",79);
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::init ()
{
  this->initSolution(this->getNoDOFs(), 1 + drc.getOrder());

  this->initSystem(Dim::opt.solver);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);
  if (mVec.size() == 1)
    drc.setMaterial(mVec.front());

  return drc.getOrder() > 0;
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::solveStep (const TimeStep& tp, bool forceNewTan)
{
  // Schedule change, reinit buffers
  if (!this->setMode(SIM::DYNAMIC))
    return false;

  if (tp.step == drc.getOrder() && drc.lCache())
    this->initLHSbuffers();

  Vector dummy;
  this->updateDirichlet(tp.time.t,&dummy);
  if (!this->assembleSystem(tp.time, solution, newTangent || forceNewTan))
    return false;

  if (!this->solveSystem(solution.front()))
    return false;

  this->printSolutionSummary(solution.front(), 0, "concentration");

  newTangent = tp.step < drc.getOrder() || !drc.lCache();

  return true;
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::advanceStep (TimeStep&)
{
  drc.advanceStep();
  this->pushSolution();

  return true;
}


template<class Dim>
bool SIMDarcyAdvection<Dim>::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size())
    return false;

  drc.setMaterial(mVec[propInd]);

  return true;
}


template<class Dim>
int SolverConfigurator<SIMDarcyAdvection<Dim>>::
setup (SIMDarcyAdvection<Dim>& darcy,
       const typename SIMDarcyAdvection<Dim>::SetupProps& props,
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


//! \brief Instantiation macro.
#define INSTANTIATE(T) \
  template class SIMDarcyAdvection<T>; \
  template struct SolverConfigurator<SIMDarcyAdvection<T>>;

INSTANTIATE(SIM1D)
INSTANTIATE(SIM2D)
INSTANTIATE(SIM3D)
