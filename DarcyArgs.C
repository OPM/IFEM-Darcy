// $Id$
//==============================================================================
//!
//! \file DarcyArgs.C
//!
//! \date Mar 7 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Darcy application.
//!
//==============================================================================

#include "DarcyArgs.h"
#include "Utilities.h"
#include "tinyxml.h"


bool DarcyArgs::parse(const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"darcy"))
    utl::getAttribute(elem,"adaptive",adap);

  return SIM::AppXMLInputBase::parse(elem);
}
