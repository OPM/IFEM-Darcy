// $Id$
//==============================================================================
//!
//! \file DarcyMaterial.h
//!
//! \date Oct 24 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Material implementations for Darcy flow problems.
//!
//==============================================================================

#ifndef _DARCY_MATERIAL_H_
#define _DARCY_MATERIAL_H_

#include <memory>

class RealFunc;
namespace tinyxml2 { class XMLElement; }
class VecFunc;
class Vec3;

/*!
  \brief Class representing the integrand of the Darcy problem.
*/

class DarcyMaterial
{
public:
  const double rhow = 1.0; //!< Density of fluid

  //! \brief Check if a tag is for material data.
  static bool handlesTag(const char* name);

  //! \brief Parse an XML block.
  bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Sets the permeability function.
  void setPermValues(VecFunc* perm);
  //! \brief Sets the permeability scalar field function.
  void setPermField(RealFunc* perm);
  //! \brief Returns the permeability at a given point.
  Vec3 getPermeability(const Vec3& X) const;

  //! \brief Returns the dispersivity at a given point.
  double getDispersivity(const Vec3& X) const;

  //! \brief Defines a scalar dispersivity function.
  void setDispersivity(RealFunc* f);

  //! \brief Returns the porosity at a given point.
  double getPorosity(const Vec3& X) const;

  //! \brief Defines a scalar porosity function.
  void setPorosity(RealFunc* f);

  //! \brief Define fluid viscosity used to calculate mobilities.
  void setViscosity(double nu) { viscosity = nu; }

  //! \brief Obtain fluid viscosity.
  double getViscosity() const { return viscosity; }

  //! \brief Helper for CoSTA.
  void setParam(const std::string& name, double value);

protected:
  std::unique_ptr<VecFunc>  permvalues;   //!< Permeability function
  std::unique_ptr<RealFunc> permeability; //!< Permeability field function
  std::unique_ptr<RealFunc> porosity;     //!< Porosity function
  std::unique_ptr<RealFunc> dispersivity; //!< Dispersivity function
  double viscosity{}; //!< Fluid viscosity
};

#endif
