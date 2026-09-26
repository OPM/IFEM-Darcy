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
#include <cstring>


DarcyMaterial::DarcyMaterial (const tinyxml2::XMLElement* elem)
{
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    this->parse(child);
}


DarcyMaterial::DarcyMaterial (DarcyMaterial&& tmp)
{
  if (!permmatrix.get())
    permmatrix = std::move(tmp.permmatrix);

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

  if (viscosity == 1.0 && tmp.viscosity > 0.0)
    viscosity = tmp.viscosity;
  else if (tmp.viscosity < 0.0)
    viscosity = tmp.viscosity; // invalid properties
}


DarcyMaterial::DarcyMaterial () = default;


DarcyMaterial::~DarcyMaterial () = default;


/*!
  \return \e false if \a elem is not a material tag, otherwise \a true
  (also if the material data is invalid, which is flagged by setting
  \ref viscosity to -1.0)
*/

bool DarcyMaterial::parse (const tinyxml2::XMLElement* elem)
{
  std::string type;
  utl::getAttribute(elem,"type",type);

  const char* value = nullptr;
  if (utl::getValue(elem,"permvalues"))
  {
    std::cerr <<" *** DarcyMaterial::parse(): The <permvalues> tag is not"
              <<" supported, use <permeability type=\"diag\"> instead."
              << std::endl;
    viscosity = -1.0;
  }
  else if ((value = utl::getValue(elem,"permeability")))
  {
    IFEM::cout <<"\t\tPermeability";
    if (type == "diag" || type == "diagonal" || type == "vector")
      permvalues.reset(utl::parseVecFunc(value));
    else if (type != "matrix")
      permeability.reset(utl::parseRealFunc(value,type));
    else
    {
      std::vector<double> values;
      char* strtmp = strdup(value);
      for (char* s = strtok(strtmp," \t\r\n"); s; s = strtok(nullptr," \t\r\n"))
        values.push_back(atof(s));
      free(strtmp);
      size_t ndim = static_cast<size_t>(sqrt(values.size()));
      if (ndim < 2 || ndim > 3)
      {
        IFEM::cout << std::endl;
        std::cerr <<" *** DarcyMaterial::parse(): Invalid matrix dimension ("
                  << ndim <<"), "<< values.size() <<" values specified."
                  << std::endl;
        viscosity = -1.0;
      }
      else
      {
        Matrix Kmat(ndim,ndim);
        Kmat.fill(values.data());
        permmatrix.reset(new Matrix(Kmat,true));
        IFEM::cout <<" matrix:"<< *permmatrix;
      }
      return true;
    }
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
  else if (!(value = utl::getValue(elem,"viscosity")))
    return false; // not a material tag
  else if (double mu = atof(value); mu > 0.0 && viscosity > 0.0)
    IFEM::cout <<"\t\tFluid viscosity: "<< (viscosity = mu) << std::endl;

  return true;
}


Vec3 DarcyMaterial::getPermeability (const Vec3& X) const
{
  Vec3 result;
  if (permvalues.get())
    result = (*permvalues)(X);
  else if (permeability.get())
    result = (*permeability)(X);
  else
    result = 1.0;

  return result;
}


bool DarcyMaterial::getPermeability (Matrix* K) const
{
  if (!permmatrix.get())
    return false;
  else if (K)
    *K = *permmatrix;

  return true;
}


double DarcyMaterial::getPorosity (const Vec3& X) const
{
  return porosity.get() ? (*porosity)(X) : 0.0;
}


double DarcyMaterial::getDispersivity (const Vec3& X) const
{
  return dispersivity.get() ? (*dispersivity)(X) : 0.0;
}


double DarcyMaterial::getDensity (double c) const
{
  double rho = density.get() ? (*density)(c) : 1.0;
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
