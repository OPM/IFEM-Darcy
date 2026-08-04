// $Id$
//==============================================================================
//!
//! \file SIMDarcySchedule.h
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for scheduled Darcy advection problems.
//!
//==============================================================================

#ifndef _SIM_DARCY_SCHEDULE_H_
#define _SIM_DARCY_SCHEDULE_H_

#include "SIMadmin.h"
#include "SIMDarcy.h"
#include "SIMDarcyAdvection.h"
#include "SIMCoupled.h"
#include <set>


/*!
  \brief Driver class for analysis of scheduled Darcy advection problems.
*/

template<class Dim>
class SIMDarcySchedule : public SIMCoupled<SIMDarcy<Dim>, SIMDarcyAdvection<Dim>>,
                         public SIMadmin
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMDarcySchedule(SIMDarcy<Dim>& s1, SIMDarcyAdvection<Dim>& s2) :
    SIMCoupled<SIMDarcy<Dim>,SIMDarcyAdvection<Dim>>(s1,s2),
    schit(schedule.end()) {}

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool = true) override;

  //! \brief Sets up field dependencies.
  void setupDependencies() override;

  //! \brief Parses an XML input element.
  bool parse(const tinyxml2::XMLElement* elem) override;

protected:
  std::set<double> schedule; //!< Scheduled pressure changes
  std::set<double>::const_iterator schit; //!< Iterator for current schedule entry
};

#endif
