// $Id$
//==============================================================================
//!
//! \file DarcyMaterial.C
//!
//! \date Oct 24 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Material implementations for Darcy flow problems.
//!
//==============================================================================

#include "DarcyMaterial.h"
#include "Functions.h"
#include "IFEM.h"
#include "Utilities.h"
#include "Vec3.h"
#include "tinyxml2.h"


DarcyMaterial::DarcyMaterial (const tinyxml2::XMLElement* elem)
{
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    this->parse(child);
}


DarcyMaterial::DarcyMaterial (DarcyMaterial&& tmp)
{
  if (!permvalues.get())
    permvalues = std::move(tmp.permvalues);

  if (!permeability.get())
    permeability = std::move(tmp.permeability);

  if (!porosity.get())
    porosity = std::move(tmp.porosity);

  if (!dispersivity.get())
    dispersivity = std::move(tmp.dispersivity);

  if (!density.get())
    density = std::move(tmp.density);

  if (viscosity == 0.0)
    viscosity = tmp.viscosity;
}


DarcyMaterial::DarcyMaterial () = default;


DarcyMaterial::~DarcyMaterial () = default;


bool DarcyMaterial::parse (const tinyxml2::XMLElement* elem)
{
  std::string type;
  utl::getAttribute(elem,"type",type);

  const char* value = nullptr;
  if ((value = utl::getValue(elem,"permvalues")))
  {
    IFEM::cout <<"\t\tPermeability";
    permvalues.reset(utl::parseVecFunc(value));
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"permeability")))
  {
    IFEM::cout <<"\t\tPermeability";
    permeability.reset(utl::parseRealFunc(value,type));
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"porosity")))
  {
    IFEM::cout <<"\t\tPorosity";
    porosity.reset(utl::parseRealFunc(value,type));
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"dispersivity")))
  {
    IFEM::cout <<"\t\tDispersivity";
    dispersivity.reset(utl::parseRealFunc(value,type));
    IFEM::cout << std::endl;
  }
  else if ((value = utl::getValue(elem,"density")))
  {
    IFEM::cout <<"\t\tFluid density: ";
    density.reset(utl::parseTimeFunc(value,type));
  }
  else if ((value = utl::getValue(elem,"viscosity")))
    IFEM::cout <<"\t\tFluid viscosity: "
               << (viscosity = atof(value)) << std::endl;
  else
    return false;

  return true;
}


Vec3 DarcyMaterial::getPermeability (const Vec3& X) const
{
  Vec3 result;
  if (permvalues)
    result = (*permvalues)(X);
  else if (permeability)
    result = (*permeability)(X);

  return result;
}


double DarcyMaterial::getPorosity (const Vec3& X) const
{
  return porosity ? (*porosity)(X) : 0.0;
}


double DarcyMaterial::getDispersivity (const Vec3& X) const
{
  return dispersivity ? (*dispersivity)(X) : 0.0;
}


double DarcyMaterial::getDensity (double c) const
{
  double rho = density ? (*density)(c) : 1.0;
  if (rho > 1.0e-16) return rho;

  std::cerr <<" *** DarcyMaterial::getDensity(): Non-positive fluid density ("
            << rho <<")"<< std::endl;
  return -1.0;
}


void DarcyMaterial::setParam (const std::string& name, double value)
{
  if (permvalues.get())
    permvalues->setParam(name,value);

  if (permeability.get())
    permeability->setParam(name,value);

  if (porosity.get())
    porosity->setParam(name,value);

  if (dispersivity.get())
    dispersivity->setParam(name,value);
}
