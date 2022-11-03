// $Id$
//==============================================================================
//!
//! \file SIMDarcySchedule.C
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of scheduled Darcy advection.
//!
//==============================================================================

#include "SIMDarcySchedule.h"

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"
#include <tinyxml.h>


template<class Dim>
SIMDarcySchedule<Dim>::SIMDarcySchedule(SIMDarcy<Dim>& dcySim,
                                        SIMDarcyAdvection<Dim>& advSim) :
  Base(dcySim, advSim)
{
}


template<class Dim>
bool SIMDarcySchedule<Dim>::solveStep (TimeStep& tp, bool)
{
  if (schedule.empty()) {
    std::cerr << "*** No schedule configured. Bailing\n";
    return false;
  }

  pressureSolved = false;

  if (tp.step == 1) {
    if (!this->S1.solveStep(tp))
      return false;
    pressureSolved = true;
  } else if (currSchedule < schedule.size() &&
      tp.time.t >= schedule[currSchedule]) {
    IFEM::cout << "\n* Scheduled pressure change at time "
               << schedule[currSchedule]
                  <<", calculating new pressure matrix." << std::endl;
    if (!this->S1.solveStep(tp))
      return false;

    ++currSchedule;
    this->S2.newTangent = true;
    pressureSolved = true;
  } else
    this->S1.printStep(tp.step, tp.time);

  return this->S2.solveStep(tp);
}

template<class Dim>
bool SIMDarcySchedule<Dim>::saveStep(const TimeStep& tp, int& nBlock)
{
  return this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock,pressureSolved);
}


template<class Dim>
void SIMDarcySchedule<Dim>::setupDependencies ()
{
  this->S2.registerDependency(&this->S1, "pressure", 1, this->S1.getFEModel(), 1);
}


template<class Dim>
bool SIMDarcySchedule<Dim>::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"darcy"))
    return true;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    if (!strcasecmp(child->Value(),"schedule")) {
      const TiXmlElement* sched = child->FirstChildElement("update");
      IFEM::cout << "\n\tScheduled pressure updates:";
      for (; sched; sched = sched->NextSiblingElement("update")) {
        double t;
        if (utl::getAttribute(sched,"time",t)) {
          IFEM::cout << "\n\t\tUpdate at " << t;
          schedule.push_back(t);
        }
      }
      IFEM::cout << std::endl;
    }
  }

  return true;
}


template class SIMDarcySchedule<SIM1D>;
template class SIMDarcySchedule<SIM2D>;
template class SIMDarcySchedule<SIM3D>;
