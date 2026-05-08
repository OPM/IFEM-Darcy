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


/*!
  \brief Driver class for analysis of scheduled Darcy advection problems.
*/

template<class Dim>
class SIMDarcySchedule : public SIMCoupled<SIMDarcy<Dim>, SIMDarcyAdvection<Dim>>,
                         public SIMadmin
{
  //! Convenience type alias
  using Base = SIMCoupled<SIMDarcy<Dim>,SIMDarcyAdvection<Dim>>;

public:
  //! \brief Default constructor.
  SIMDarcySchedule(SIMDarcy<Dim>& dcySim, SIMDarcyAdvection<Dim>& advSim);

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool = true) override;

  //! \brief Sets up field dependencies.
  void setupDependencies() override;

  //! \brief Parse an XML input element.
  bool parse(const tinyxml2::XMLElement* elem) override;

protected:
  size_t currSchedule = 0; //!< Index for current schedule entry
  std::vector<double> schedule; //!< Scheduled pressure changes
};

#endif
