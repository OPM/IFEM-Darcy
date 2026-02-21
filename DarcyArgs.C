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
#include "ASMmxBase.h"
#include "Utilities.h"

#include <string>
#include <strings.h>
#include "tinyxml2.h"


bool DarcyArgs::parseArg (const char* argv)
{
  if (argv[0] != '-')
    return false;
  else if (strcasecmp(argv, "-tracer") == 0)
    tracer = true;
  else if (strcmp(argv, "-mixed1") == 0)
    mixed = 1;
  else if (strcmp(argv, "-mixed") == 0)
    mixed = 2;
  else if (strcmp(argv, "-Mixed1") == 0)
    mixed = 11;
  else if (strcmp(argv, "-Mixed") == 0)
    mixed = 12;
  else if ((timeMethod = TimeIntegration::get(argv+1)) > TimeIntegration::NONE)
    ;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool DarcyArgs::parseArgComplex (int argc, char** argv, int& i)
{
  if (strcasecmp(argv[i],"-adap"))
    return false;

  adap = true;
  int j = i+1;
  if (j >= argc || argv[j][0] == '-')
    adNorm = DCY::PRESSURE_H1;
  else if (strcasecmp(argv[j], "pressure") == 0)
    adNorm = DCY::PRESSURE_H1, ++i;
  else if (strcasecmp(argv[j], "recovery_press") == 0)
    adNorm = DCY::RECOVERY_PRESSURE, ++i;
  else if (strcasecmp(argv[j], "concentration") == 0)
    adNorm = DCY::CONCENTRATION_H1, ++i;
  else if (strcasecmp(argv[j], "recovery_conc") == 0)
    adNorm = DCY::RECOVERY_CONCENTRATION, ++i;
  else if (strcasecmp(argv[j], "total") == 0)
    adNorm = DCY::TOTAL_H1, ++i;
  else if (strcasecmp(argv[j], "recovery") == 0)
    adNorm = DCY::RECOVERY, ++i;
  else
    adNorm = DCY::PRESSURE_H1;

  return true;
}


bool DarcyArgs::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    std::string type;
    if (utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
  }
  else if (!strcasecmp(elem->Value(),"darcy")) {
    utl::getAttribute(elem,"tracer",tracer);
    ASMmxBase::Type = ASMmxBase::NONE;
    if (const char* ad = elem->Attribute("adap"); ad) {
      if (strcasecmp(ad,"pressure") == 0)
        adNorm = DCY::PRESSURE_H1;
      else if (strcasecmp(ad,"recovery_press") == 0)
        adNorm = DCY::RECOVERY_PRESSURE;
      else if (strcasecmp(ad,"concentration") == 0)
        adNorm = DCY::CONCENTRATION_H1;
      else if (strcasecmp(ad,"recovery_conc") == 0)
        adNorm = DCY::RECOVERY_CONCENTRATION;
      else if (strcasecmp(ad,"total") == 0)
        adNorm = DCY::TOTAL_H1;
      else if (strcasecmp(ad,"recovery") == 0)
        adNorm = DCY::RECOVERY;
    }

    if (elem->FirstChildElement("schedule"))
      scheduled = true;
  }

  return this->SIMargsBase::parse(elem);
}
