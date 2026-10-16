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
  //! \brief Enumeration of adaptation norms for Darcy problems.
  enum AdaptationNorm {
    NO_ADAP,     //!< No mesh adaptation specified
    TOTAL_H1,    //!< Based on total H1 error
    PRESSURE_H1, //!< Based on H1 error of pressure (or L2 for div-compatible)
    CONCENTR_H1, //!< Based on H1 error of concentration
    FLUX_H1,     //!< Based on H1 error of Darcy flux (div-compatible only)
    RECOVERY,          //!< Based on total recovered H1 error
    RECOVERY_PRESSURE, //!< Based on H1 error of recovered pressure
    RECOVERY_CONCENTR  //!< Based on H1 error of recovered concentration
  };
}

#endif
