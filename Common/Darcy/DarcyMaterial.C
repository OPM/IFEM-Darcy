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

#include <tinyxml2.h>


bool DarcyMaterial::handlesTag(const char* name)
{
  return !strcasecmp(name, "dispersivity") ||
         !strcasecmp(name, "permvalues") ||
         !strcasecmp(name, "permeability") ||
         !strcasecmp(name, "porosity") ||
         !strcasecmp(name, "viscosity");
}


bool DarcyMaterial::parse (const tinyxml2::XMLElement* elem)
{
  auto doParse = [this](const tinyxml2::XMLElement* child)
  {
    const char* value = nullptr;
    if ((value = utl::getValue(child,"permvalues"))) {
      IFEM::cout <<"\t\tPermeability";
      this->setPermValues(utl::parseVecFunc(value));
      IFEM::cout << std::endl;
    } else if ((value = utl::getValue(child,"permeability"))) {
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\t\tPermeability";
      this->setPermField(utl::parseRealFunc(value,type));
      IFEM::cout << std::endl;
    } else if ((value = utl::getValue(child,"porosity"))) {
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\t\tPorosity";
      this->setPorosity(utl::parseRealFunc(value,type));
      IFEM::cout << std::endl;
    } else if ((value = utl::getValue(child,"dispersivity"))) {
      std::string type;
      utl::getAttribute(child,"type",type);
      IFEM::cout <<"\t\tDispersivity";
      this->setDispersivity(utl::parseRealFunc(value,"expression"));
      IFEM::cout << std::endl;
    } else if ((value = utl::getValue(child,"viscosity"))) {
      double viscosity = atof(value);
      IFEM::cout << "\t\tFluid viscosity: " << viscosity << std::endl;
      this->setViscosity(viscosity);
    } else
      return false;

    return true;
  };

  if (strcasecmp(elem->Value(), "materialdata"))
    return doParse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!doParse(child))
      return false;

  return true;
}


void DarcyMaterial::setPermValues (VecFunc* perm)
{
  permvalues.reset(perm);
}


void DarcyMaterial::setPermField (RealFunc* perm)
{
  permeability.reset(perm);
}


void DarcyMaterial::setDispersivity (RealFunc* f)
{
  dispersivity.reset(f);
}


void DarcyMaterial::setPorosity (RealFunc* f)
{
  porosity.reset(f);
}


double DarcyMaterial::getDispersivity(const Vec3& X) const
{
  return dispersivity ? (*dispersivity)(X) : 0.0;
}


double DarcyMaterial::getPorosity (const Vec3& X) const
{
  return porosity ? (*porosity)(X) : 0.0;
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


void DarcyMaterial::setParam (const std::string& name, double value)
{
  if (permvalues.get())
    permvalues->setParam(name,value);

  if (porosity.get())
    porosity->setParam(name,value);

  if (dispersivity.get())
    dispersivity->setParam(name,value);

  if (permeability.get())
    permeability->setParam(name,value);
}
