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
#include "tinyxml.h"


bool DarcyArgs::parseArg (const char* argv)
{
  TimeIntegration::Method tmp;
  if (argv[0] != '-')
    return false;
  else if (strcasecmp(argv, "-twofield") == 0)
    twofield = true;
  else if ((tmp = TimeIntegration::get(argv+1)) > TimeIntegration::NONE)
    timeMethod = tmp;
  else
    return this->SIMargsBase::parseArg(argv);

  return true;
}


bool DarcyArgs::parseArgComplex (int argc, char** argv, int& i)
{
  if (!strcmp(argv[i],"-adap")) {
    adNorm = DCY::PRESSURE_H1;
    adap = true;
    if (i+1 < argc && argv[i+1][0] != '-' &&
        strcasecmp(argv[i+1], "pressure") == 0)
      adNorm = DCY::PRESSURE_H1, ++i;
    else if (i+1 < argc && argv[i+1][0] != '-' &&
             strcasecmp(argv[i+1], "recovery_press") == 0)
      adNorm = DCY::RECOVERY_PRESSURE, ++i;
    else if (i+1 < argc && argv[i+1][0] != '-' &&
             strcasecmp(argv[i+1], "concentration") == 0)
      adNorm = DCY::CONCENTRATION_H1, ++i;
    else if (i+1 < argc && argv[i+1][0] != '-' &&
             strcasecmp(argv[i+1], "recovery_conc") == 0)
      adNorm = DCY::RECOVERY_CONCENTRATION, ++i;
    else if (i+1 < argc && argv[i+1][0] != '-' &&
             strcasecmp(argv[i+1], "total") == 0)
      adNorm = DCY::TOTAL_H1, ++i;
    else if (i+1 < argc && argv[i+1][0] != '-' &&
             strcasecmp(argv[i+1], "recovery") == 0)
      adNorm = DCY::RECOVERY, ++i;
  } else
    return false;

  return true;
}


bool DarcyArgs::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"timestepping")) {
    std::string type;
    if (utl::getAttribute(elem,"type",type))
      timeMethod = TimeIntegration::get(type);
  }
  if (!strcasecmp(elem->Value(),"darcy")) {
    utl::getAttribute(elem,"twofield",twofield);
    ASMmxBase::Type = ASMmxBase::NONE;
    const char* ad = elem->Attribute("adap");
    if (ad) {
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
  }

  return this->SIMargsBase::parse(elem);
}
