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
class ScalarFunc;
class VecFunc;
class Vec3;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class representing a physical property item for Darcy flow problems.
*/

class DarcyMaterial
{
public:
  //! \brief Constructor parsing properties from an XML-element.
  explicit DarcyMaterial(const tinyxml2::XMLElement* elem);
  //! \brief Move constructor.
  DarcyMaterial(DarcyMaterial&& tmp);
  //! \brief Default constructor.
  DarcyMaterial();
  //! \brief Default destructor.
  virtual ~DarcyMaterial();

  //! \brief Parses an XML-element.
  bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Returns the permeability at a given point.
  Vec3 getPermeability(const Vec3& X) const;

  //! \brief Returns the porosity at a given point.
  double getPorosity(const Vec3& X) const;

  //! \brief Returns the dispersivity at a given point.
  double getDispersivity(const Vec3& X) const;

  //! \brief Returns the fluid density.
  double getDensity(double c) const;

  //! \brief Returns the fluid viscosity.
  double getViscosity() const { return viscosity; }

  //! \brief Helper for CoSTA.
  void setParam(const std::string& name, double value);

private:
  std::unique_ptr<VecFunc>    permvalues;   //!< Permeability function
  std::unique_ptr<RealFunc>   permeability; //!< Permeability field function
  std::unique_ptr<RealFunc>   porosity;     //!< Porosity function
  std::unique_ptr<RealFunc>   dispersivity; //!< Dispersivity function
  std::unique_ptr<ScalarFunc> density;      //!< Fluid density function

  double viscosity = 0.0; //!< Fluid viscosity
};

#endif
