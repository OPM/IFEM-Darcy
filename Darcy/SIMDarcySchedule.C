// $Id$
//==============================================================================
//!
//! \file SIMDarcySchedule.C
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for scheduled Darcy advection problems.
//!
//==============================================================================

#include "SIMDarcySchedule.h"

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"
#include <tinyxml2.h>


template<class Dim>
bool SIMDarcySchedule<Dim>::solveStep (TimeStep& tp, bool)
{
  bool newTangent = false;
  bool pressureOk = true;
  if (tp.step == 1)
    pressureOk = this->S1.solveStep(tp);
  else if (schit != schedule.end() && tp.time.t >= *schit)
  {
    IFEM::cout <<"\n * Scheduled pressure change at time "<< *schit
               <<", calculating new pressure matrix."<< std::endl;
    newTangent = pressureOk = this->S1.solveStep(tp);
    while (schit != schedule.end() && *schit <= tp.time.t) ++schit;
  }
  else
    this->S1.keepStep(tp);

  return pressureOk && this->S2.solveStep(tp,newTangent);
}


template<class Dim>
void SIMDarcySchedule<Dim>::setupDependencies ()
{
  this->S1.registerDependency(&this->S2,"tracer"  ,1,this->S2.getFEModel(),1);
  this->S2.registerDependency(&this->S1,"pressure",1,this->S1.getFEModel(),1);
}


template<class Dim>
bool SIMDarcySchedule<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"darcy"))
    return true;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"schedule"))
    {
      const tinyxml2::XMLElement* sched = child->FirstChildElement("update");
      for (; sched; sched = sched->NextSiblingElement("update"))
        if (double t = 0.0; utl::getAttribute(sched,"time",t))
          schedule.insert(t);
        else if (utl::getAttribute(sched,"start",t))
          if (double t1 = 0.0; utl::getAttribute(sched,"end",t1) && t1 > t)
            if (double dt = 0.0; utl::getAttribute(sched,"dt",dt) && dt > 0.0)
              for (; t <= t1; t += dt)
                schedule.insert(t);
      schit = schedule.begin();
      IFEM::cout <<"\tScheduled pressure updates:";
      size_t i = 0;
      for (double t : schedule)
        IFEM::cout << (schedule.size() > 20 && (++i)%20 == 1 ? "\n\t" : " ") << t;
      IFEM::cout << std::endl;
    }

  return true;
}


template class SIMDarcySchedule<SIM1D>;
template class SIMDarcySchedule<SIM2D>;
template class SIMDarcySchedule<SIM3D>;
