// $Id$
//==============================================================================
//!
//! \file DarcyArgs.C
//!
//! \date Sep 15 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Darcy application.
//!
//==============================================================================

#include "DarcyArgs.h"
#include "Utilities.h"
#include "tinyxml.h"


bool DarcyArgs::parseArg (const char* argv)
{
  TimeIntegration::Method tmp;
  if (argv[0] != '-')
    return false;
  else if ((tmp = TimeIntegration::get(argv+1)) > TimeIntegration::NONE)
    timeMethod = tmp;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool DarcyArgs::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    std::string type;
    if (utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
  }

  return this->SIMargsBase::parse(elem);
}
