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

#include "SIMargsBase.h"
#include "TimeIntUtils.h"


class TiXmlElement;


/*!
  \brief Class holding application parameters.
*/

class DarcyArgs : public SIMargsBase
{
public:
  TimeIntegration::Method timeMethod = TimeIntegration::NONE; //!< Time integration method

  //! \brief Default constructor.
  DarcyArgs() : SIMargsBase("darcy") {}

  //! \brief Parses a command-line argument.
  bool parseArg(const char* argv) override;

protected:
  //! \brief Parse an element from the input file
  bool parse(const TiXmlElement* elem) override;
};

#endif
