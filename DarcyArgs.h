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
  TimeIntegration::Method timeMethod = TimeIntegration::NONE; //!< Time integration method
  bool tracer = false; //!< Include a tracer field
  bool scheduled = false; //!< If true, update pressure on a schedule, separate ad step.
  DCY::AdaptationNorm adNorm = DCY::NO_ADAP; //!< Norm to adapt based on

  //! \brief Default constructor.
  DarcyArgs() : SIMargsBase("darcy") {}

  //! \brief Parses a command-line argument.
  bool parseArg(const char* argv) override;

  //! \brief Parses a command-line argument with parameters.
  bool parseArgComplex(int argc, char **argv, int &i);

protected:
  //! \brief Parse an element from the input file
  bool parse(const tinyxml2::XMLElement* elem) override;
};

#endif
