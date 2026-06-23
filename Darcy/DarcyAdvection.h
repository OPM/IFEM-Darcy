// $Id$
//==============================================================================
//!
//! \file DarcyAdvection.h
//!
//! \date Aug 26 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy advection problems.
//!
//==============================================================================

#ifndef _DARCY_ADVECTION_H_
#define _DARCY_ADVECTION_H_

#include "DarcyBase.h"

#include <memory>

class Field;
class RealFunc;
class VecFunc;


/*!
  \brief Class representing the integrand of the Darcy transport problem.
*/

class DarcyAdvection : public DarcyBase
{
public:
  //! \brief The constructor initializes all pointers to zero.
  explicit DarcyAdvection(unsigned short int n, int torder = 0);
  //! \brief The destructor deletes the \ref bodyforce.
  virtual ~DarcyAdvection();

  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  using DarcyBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  using DarcyBase::evalSol2;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] eV Element vectors at current point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? nsd+1  : 1; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Index for field
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t i, const char* prefix) const override;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  NormBase* getNormIntegrand(AnaSol*) const override;

  //! \brief Set a field from a dependency.
  void setNamedField(const std::string& name, Field* field) override;

protected:
  VecFunc* bodyforce; //!< Body force function

  std::unique_ptr<RealFunc> source; //!< Tracer source function
  std::unique_ptr<Field>    pField; //!< Pressure field
  Vector                    pVec;   //!< Pressure values

  friend class DarcyAdvectionNorm;
};


/*!
  \brief Class representing the integrand of Darcy advection energy norms.
*/

class DarcyAdvectionNorm : public NormBase
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit DarcyAdvectionNorm(DarcyAdvection& p) : NormBase(p) {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t, size_t j, const char* prefix) const override;
};

#endif
