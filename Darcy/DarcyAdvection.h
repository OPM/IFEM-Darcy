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

#include "DarcyMaterial.h"

#include "BDF.h"
#include "IntegrandBase.h"

#include <memory>

class Field;
class RealFunc;
class SIMbase;


/*!
  \brief Class representing the integrand of the Darcy transport problem.
*/

class DarcyAdvection : public IntegrandBase
{
public:
  //! \brief The constructor initializes all pointers to zero.
  DarcyAdvection(unsigned short int n, int torder = 0);
  //! \brief Empty destructor.
  virtual ~DarcyAdvection();

  //! \brief Assigns the owner simulator (used by parse()).
  void setOwnerSim(SIMbase* sim) { ownerSim = sim; }

  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Element number
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                  bool neumann) const override;

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] XC Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points on this element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const FiniteElement& fe, const Vec3& XC, size_t nPt,
                   LocalIntegral& elmInt) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  using IntegrandBase::evalSol2;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] eV Element vectors at current point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] iGP Global integration point counter of first point in element
  bool finalizeElement(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, size_t iGP) override;

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

  //! \brief Returns order of time integration.
  int getOrder() const { return bdf.getActualOrder(); }

  //! \brief Update time stepping scheme (BE -> BDF2 transition).
  void advanceStep()
  {
    bdf.advanceStep();
  }

  //! \brief Set a field from a dependency.
  void setNamedField(const std::string& name, Field* field) override;

  //! \brief Initializes and toggles the use of left-hand-side matrix buffers.
  //! \param[in] nEl Number of elements in the model/toggle.
  //! - If larger than 1, element matrix buffers are allocated to given size.
  //! - If equal to 1, element matrices are recomputed.
  //! - If equal to 0, reuse cached element matrices.
  void initLHSbuffers(size_t nEl) override;

  //! \brief Enable/disable caching of element matrices.
  void lCache(bool enable) { calcMats = enable ? 1 : -1; }
  //! \brief Returns whether or not caching of element matrices is enabled.
  bool lCache() const { return calcMats >= 0; }

  //! \brief Set material parameters.
  void setMaterial(DarcyMaterial& mat1) { mat = &mat1; }

protected:
  SIMbase* ownerSim; //!< The simulator that owns this integrand

  const DarcyMaterial* mat; //!< Material to use

  std::unique_ptr<RealFunc> source; //!< Tracer source function
  std::unique_ptr<Field>    pField; //!< Pressure field
  Vector                    pVec;   //!< Pressure values

  TimeIntegration::BDF bdf; //!< BDF time stepping helper

  Matrices myKmats; //!< Cached element matrices

  //! Flag for calculation/caching of element matrices.
  //! - < 0 : Always recalculate the element matrices
  //! - = 0 : Use cached element matrices
  //! - > 1 : Calculate new element matrices
  char calcMats = -1;

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
