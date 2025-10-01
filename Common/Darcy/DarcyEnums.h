// $Id$
//==============================================================================
//!
//! \file DarcyEnums.h
//!
//! \date Oct 27 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Some enumerations used in Darcy problems.
//!
//==============================================================================
#ifndef _DARCY_ENUMS_H_
#define _DARCY_ENUMS_H_


namespace DCY //!< Darcy scope
{
  //! \brief Enumeration of Darcy adaptation norms.
  enum AdaptationNorm {
    NO_ADAP, //!< No adaptation
    TOTAL_H1, //!< Adapt based on total H1 error
    PRESSURE_H1, //!< Adapt based on H1 error of pressure
    CONCENTRATION_H1, //!< Adapt based on H1 error of concentration
    RECOVERY, //!< Adapt based on total recovered H1 error
    RECOVERY_PRESSURE, //!< Adapt based on H1 error of recovered pressure
    RECOVERY_CONCENTRATION, //!< Adapt based on H1 error of recovered concentration
  };
}

#endif
