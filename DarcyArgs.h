// $Id$
//==============================================================================
//!
//! \file DarcyArgs.h
//!
//! \date Jan 24 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Darcy application.
//!
//==============================================================================
#ifndef DARCY_ARGS_H
#define DARCY_ARGS_H

#include "AppCommon.h"


/*! \brief Struct holding application parameters.
 */

class DarcyArgs : public SIM::AppXMLInputBase
{
public:
  bool adap = false; //!< True to run an adaptive simulator

protected:
  //! \brief Parse an element from the input file
  bool parse(const TiXmlElement* elem) override;
};

#endif
