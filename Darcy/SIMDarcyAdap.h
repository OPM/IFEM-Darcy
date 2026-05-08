// $Id$
//==============================================================================
//!
//! \file SIMDarcyAdap.h
//!
//! \date Oct 27 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for adaptive NURBS-based FE analysis of Darcy problems.
//!
//==============================================================================
#ifndef _SIM_DARCY_ADAP_H
#define _SIM_DARCY_ADAP_H

#include "DarcyEnums.h"
#include "SIMSolverAdap.h"


/*!
  \brief Adaptive simulator for Darcy handling norm indices.
*/

template<class T1>
class SIMDarcyAdap : public SIMSolverAdapImpl<T1, AdaptiveISolver<T1>>
{
public:
  //! \brief Default constructor.
  //! \param s1 The darcy simulator to run for.
  explicit SIMDarcyAdap(T1& s1) : SIMSolverAdapImpl<T1, AdaptiveISolver<T1>>(s1)
  {
    switch (s1.getAdaptiveNorm()) {
      case DCY::TOTAL_H1:
        this->aSim.setAdaptationNorm(0, DarcyNorm::TOTAL_NORM_E+1); break;
      case DCY::PRESSURE_H1:
        this->aSim.setAdaptationNorm(0, DarcyNorm::H1_E_Ph+1); break;
      case DCY::CONCENTRATION_H1:
        this->aSim.setAdaptationNorm(0, DarcyNorm::H1_E_Ch+1); break;
      case DCY::RECOVERY:
        this->aSim.setAdaptationNorm(1, DarcyNorm::TOTAL_NORM_REC+1); break;
      case DCY::RECOVERY_PRESSURE:
        this->aSim.setAdaptationNorm(1, DarcyNorm::H1_Pr_Ph+1); break;
      case DCY::RECOVERY_CONCENTRATION:
        this->aSim.setAdaptationNorm(1, DarcyNorm::H1_Cr_Ch+1); break;
      default: break;
    }
  }
};

#endif
