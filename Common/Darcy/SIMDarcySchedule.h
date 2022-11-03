// $Id$
//==============================================================================
//!
//! \file SIMDarcySchedule.h
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of scheduled Darcy advection.
//!
//==============================================================================

#ifndef _SIM_DARCY_SCHEDULE_H_
#define _SIM_DARCY_SCHEDULE_H_

#include "SIMadmin.h"
#include "SIMDarcy.h"
#include "SIMDarcyAdvection.h"
#include "SIMCoupled.h"
/*!
  \brief Driver class for isogeometric FE analysis of scheduled Darcy advection problems.
*/

template<class Dim>
class SIMDarcySchedule : public SIMCoupled<SIMDarcy<Dim>, SIMDarcyAdvection<Dim>>,
                         public SIMadmin
{
  using Base = SIMCoupled<SIMDarcy<Dim>, SIMDarcyAdvection<Dim>>; //!< Convenience typedef
public:
  //! \brief Default constructor.
  SIMDarcySchedule(SIMDarcy<Dim>& dcySim, SIMDarcyAdvection<Dim>& advSim);

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp, bool=false) override;

  //! \brief Saves the converged results to VTF-file of a given time step.
  bool saveStep(const TimeStep& tp, int& nBlock) override;

  //! \brief Sets up field dependencies.
  void setupDependencies() override;

  //! \brief Parse an XML input element.
  bool parse(const TiXmlElement* elem) override;

protected:
  size_t currSchedule = 0;
  std::vector<double> schedule; //!< Scheduled pressure changes
  bool pressureSolved = false; //!< Pressure was solved for on this step
};

#endif
