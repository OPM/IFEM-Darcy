// $Id$
//==============================================================================
//!
//! \file DarcyAdvection.h
//!
//! \date Aug 26 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for tracer advection through Darcy velocity.
//!
//==============================================================================

#ifndef _DARCY_ADVECTION_H_
#define _DARCY_ADVECTION_H_

#include "DarcyMaterial.h"

#include "BDF.h"
#include "Field.h"
#include "Function.h"
#include "IntegrandBase.h"

#include "EqualOrderOperators.h"

#include <memory>

class Darcy;

/*!
  \brief Class representing the integrand of the Darcy transport problem.
*/

class DarcyAdvection : public IntegrandBase
{
  using WeakOps = EqualOrderOperators::Weak; //!< Convenience renaming

public:
  //! \brief The constructor initializes all pointers to zero.
  DarcyAdvection(unsigned short int n, const Darcy& drc, int torder = 0);
  //! \brief Empty destructor.
  virtual ~DarcyAdvection() = default;

  //! \brief Defines the tracer source function.
  void setSource(RealFunc* s) { source.reset(s); }

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iel Element number
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                  bool neumann) const override;

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points on this element
  //! \param elmInt Local integral for element
  bool initElement (const std::vector<int>& MNPC,
                    const FiniteElement& fe,
                    const Vec3& XC,
                    size_t nPt, LocalIntegral& elmInt) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time,
               const Vec3& X) const override;

  using IntegrandBase::evalSol2;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] iGP Global integration point counter of first point in element
  bool finalizeElement(LocalIntegral& elmInt,
                       const FiniteElement& fe,
                       const TimeDomain& time, size_t iGP) override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? nsd+1  : 1; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix) const override;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Returns concentration in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  double concentration(const Vectors& eV,
                       const FiniteElement& fe,
                       size_t level) const;

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
  void initLHSbuffers(size_t) override;

  //! \brief Returns whether or not to use element matrix cache.
  bool lCache() const { return useLCache; }

  //! \brief Enable/disable caching of element matrices.
  void lCache(bool enable) { useLCache = enable; }

  //! \brief Set material parameters.
  void setMaterial(DarcyMaterial& mat1) { mat = &mat1; }

protected:
  //! \brief Evaluate darcy velocity in a point.
  //! \param q Resulting darcy velocity
  //! \param fe Finite element data at current point
  //! \param X Coordinates of current point
  bool evalDarcyVel (Vector& q, const FiniteElement& fe, const Vec3& X) const;

  const DarcyMaterial* mat; //!< Material to use

  TimeIntegration::BDF bdf; //!< BDF time stepping helper

  std::unique_ptr<RealFunc> source; //!< Tracer source function

  std::unique_ptr<Field> pField; //!< Pressure field
  Vector pVec; //!< Pressure values (unused)

  const Darcy& drc; //!< Reference to darcy integrand
  bool useLCache = true; //!< True to enable caching of element matrices
  bool reuseMats = false; //!< True to reuse matrices
  Matrices myKmats; //!< Cached element matrices

  friend class DarcyAdvectionNorm;
};


/*!
  \brief Class representing the integrand of Darcy advection energy norms.
*/

class DarcyAdvectionNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] c The analytical concentration gradient (optional)
  explicit DarcyAdvectionNorm(DarcyAdvection& p);

  //! \brief Empty destructor.
  virtual ~DarcyAdvectionNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;
};

#endif
