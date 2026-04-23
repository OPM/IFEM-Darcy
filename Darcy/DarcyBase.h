// $Id$
//==============================================================================
//!
//! \file DarcyBase.h
//!
//! \date Jun 23 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for the Darcy flow integrands.
//!
//==============================================================================

#ifndef _DARCY_BASE_H_
#define _DARCY_BASE_H_

#include "HasGravityBase.h"
#include "BDF.h"

class DarcyMaterial;
class SIMbase;


/*!
  \brief Base class for the Darcy flow integrands.
*/

class DarcyBase : public HasGravityBase
{
protected:
  //! \brief The constructor is protected to allow sub-class instances only.
  DarcyBase(unsigned short int n, int torder);

public:
  //! \brief Assigns the owner simulator.
  void setOwnerSim(SIMbase* sim) { ownerSim = sim; }

  //! \brief Returns the permeability at specified point.
  Vec3 getPermeability(const Vec3& X) const;
  //! \brief Returns the dispersivity at specified point.
  double getDispersivity(const Vec3& X) const;
  //! \brief Returns the fluid viscosity.
  double getViscosity() const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  //! \brief Update time stepping scheme (BE -> BDF2 transition).
  bool advanceStep() override { return bdf.advanceStep();  }

  using HasGravityBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using HasGravityBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondence for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] XC Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points on this element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC,
                   const FiniteElement& fe, const Vec3& XC, size_t nPt,
                   LocalIntegral& elmInt) override;

  using HasGravityBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] iGP Global integration point counter of first point in element
  bool finalizeElement(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, size_t iGP) override;

  //! \brief Returns order of time integration.
  int getOrder() const { return bdf.getActualOrder(); }

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

  //! \brief Sets pointer to the material parameters object to use.
  void setMaterial(DarcyMaterial& newMat) { mat = &newMat; }

protected:
  SIMbase* ownerSim; //!< The simulator that owns this integrand

  DarcyMaterial* mat; //!< Material properties

  TimeIntegration::BDF bdf; //!< Time integration parameters

  Matrices myKmats; //!< Cached element matrices

  //! Flag for calculation/caching of element matrices.
  //! - < 0 : Always recalculate the element matrices
  //! - = 0 : Use cached element matrices
  //! - > 0 : Calculate new element matrices
  char calcMats = -1;
};

#endif
