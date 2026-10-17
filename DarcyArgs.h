// $Id$
//==============================================================================
//!
//! \file DarcyArgs.h
//!
//! \date Sep 15 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Darcy application.
//!
//==============================================================================

#ifndef _DARCY_ARGS_H
#define _DARCY_ARGS_H

#include "DarcyEnums.h"

#include "SIMargsBase.h"
#include "TimeIntUtils.h"


namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class holding application parameters.
*/

class DarcyArgs : public SIMargsBase
{
public:
  TimeIntegration::Method timeMethod; //!< Time integration method to use
  bool tracer = false; //!< Include a tracer field
  bool scheduled = false; //!< Update pressure on a schedule, separate adv. step
  bool usePiola = false; //!< Apply Piola mapping to a div-compatible flux
  bool useAL = false; //!< Use the Augmented Lagrange approach
  unsigned char mixed = 0; //!< Use a mixed formulation if &gt; 0
  DCY::AdaptationNorm adNorm = DCY::NO_ADAP; //!< Norm for mesh adaptation

  //! \brief Default constructor.
  DarcyArgs() : SIMargsBase("darcy") { timeMethod = TimeIntegration::NONE; }

  //! \brief Parses a command-line argument.
  bool parseArg(const char* argv) override;

protected:
  //! \brief Parse an element from the input file
  bool parse(const tinyxml2::XMLElement* elem) override;
};

#endif
