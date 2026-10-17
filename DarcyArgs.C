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
  // Lambda function for parsing the adaptation norm from a command-line option.
  auto parseAdapNorm = [&norm=adNorm](const char* argv)
  {
    if (strcasecmp(argv, "-pressure") == 0)
      norm = DCY::PRESSURE_H1;
    else if (strcasecmp(argv, "-recovery_press") == 0)
      norm = DCY::RECOVERY_PRESSURE;
    else if (strcasecmp(argv, "-concentration") == 0)
      norm = DCY::CONCENTR_H1;
    else if (strcasecmp(argv, "-recovery_conc") == 0)
      norm = DCY::RECOVERY_CONCENTR;
    else if (strcasecmp(argv, "-total") == 0)
      norm = DCY::TOTAL_H1;
    else if (strcasecmp(argv, "-recovery") == 0)
      norm = DCY::RECOVERY;
    else if (strcasecmp(argv, "-flux") == 0)
      norm = DCY::FLUX_H1;
    else
      norm = DCY::PRESSURE_H1;
    return true;
  };

  if (argv[0] != '-')
    return false;
  else if (strcasecmp(argv, "-tracer") == 0)
    tracer = true;
  else if (strncasecmp(argv, "-adap", 5) == 0)
    adap = parseAdapNorm(argv+5);
  else if (strcasecmp(argv, "-compatible") == 0)
    mixed = ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
  else if (strcasecmp(argv, "-piola") == 0)
    usePiola = true;
  else if (strcmp(argv, "-mixed1") == 0)
    mixed = 1;
  else if (strcmp(argv, "-mixed") == 0)
    mixed = 2;
  else if (strcmp(argv, "-Mixed1") == 0)
    mixed = 11;
  else if (strcmp(argv, "-Mixed") == 0)
    mixed = 12;
  else if (strcmp(argv, "-AL") == 0)
    useAL = true;
  else if (TimeIntegration::Method tm = TimeIntegration::get(argv+1);
           tm > TimeIntegration::NONE)
    timeMethod = tm;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool DarcyArgs::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    if (std::string type; utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
  }
  else if (!strcasecmp(elem->Value(),"darcy")) {
    utl::getAttribute(elem,"tracer",tracer);
    utl::getAttribute(elem,"piola",usePiola);
    if (std::string form; utl::getAttribute(elem,"formulation",form))
      if (form == "compatible" || form == "div-compatible")
        mixed = ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    if (const char* ad = elem->Attribute("adap"); ad) {
      if (strcasecmp(ad,"pressure") == 0)
        adNorm = DCY::PRESSURE_H1;
      else if (strcasecmp(ad,"flux") == 0)
        adNorm = DCY::FLUX_H1;
      else if (strcasecmp(ad,"recovery_press") == 0)
        adNorm = DCY::RECOVERY_PRESSURE;
      else if (strcasecmp(ad,"concentration") == 0)
        adNorm = DCY::CONCENTR_H1;
      else if (strcasecmp(ad,"recovery_conc") == 0)
        adNorm = DCY::RECOVERY_CONCENTR;
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
