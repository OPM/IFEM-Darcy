// $Id$
//==============================================================================
//!
//! \file SIMDarcyAdap.h
//!
//! \date Oct 27 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for adaptive FE analysis of Darcy problems.
//!
//==============================================================================

#ifndef _SIM_DARCY_ADAP_H
#define _SIM_DARCY_ADAP_H

#include "SIMSolverAdap.h"


/*!
  \brief Adaptive simulator for Darcy handling norm indices.
*/

template<class T1>
class SIMDarcyAdap : public SIMSolverAdapImpl<T1, AdaptiveISolver<T1>>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMDarcyAdap(T1& s1) : SIMSolverAdapImpl<T1, AdaptiveISolver<T1>>(s1)
  {
    if (unsigned int iNorm = s1.getAdaptiveNorm(); iNorm > 0)
      this->aSim.setAdaptationNorm(iNorm/100, iNorm%100);
  }

  //! \brief Returns a const reference to dummy time stepping information.
  const TimeStep& getTimePrm() const
  {
    static TimeStep dummy;
    return dummy;
  }
};

#endif
