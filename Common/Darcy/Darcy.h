// $Id$
//==============================================================================
//!
//! \file Darcy.h
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for Darcy flow problems.
//!
//==============================================================================

#ifndef _DARCY_H_
#define _DARCY_H_

#include "DarcyMaterial.h"

#include "BDF.h"
#include "EqualOrderOperators.h"
#include "IntegrandBase.h"

#include <memory>

class RealFunc;
class SIMbase;
class VecFunc;


/*!
  \brief Class representing the integrand of the Darcy problem.
*/

class Darcy : public IntegrandBase
{
protected:
  using WeakOps = EqualOrderOperators::Weak; //!< Convenience renaming

  int pp = 0; //!< Block for pressure
  int cc = 0; //!< Block for concentration
  int cp = 0; //!< Block coupling concentration to pressure

public:
  //! \brief The constructor initializes all pointers to zero.
  explicit Darcy(unsigned short int n, int torder = 0);
  //! \brief Empty destructor.
  virtual ~Darcy() {}

  //! \brief Defines the body force vector.
  void setBodyForce(VecFunc* b) { bodyforce = b; }
  //! \brief Returns the body force vector at a given point.
  Vec3 getBodyForce(const Vec3& X) const;

  //! \brief Defines the source function.
  void setSource(RealFunc* s) { source = s; }

  //! \brief Defines the concentration source function.
  virtual void setCSource(RealFunc*) {}

  //! \brief Defines a scalar flux function.
  void setFlux(RealFunc* f) { flux = f; }
  //! \brief Defines a vectorial flux function.
  void setFlux(VecFunc* f) { vflux = f; }

  //! \brief Evaluates the boundary fluid flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& normal) const;
  //! \brief Evaluates the potential source (if any) at specified point.
  double getPotential(const Vec3& X) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  //! \brief Update time stepping scheme (BE -> BDF2 transition).
  void advanceStep()
  {
    bdf.advanceStep();
  }

  //! \brief Defines the global integral for calculating reaction forces only.
  void setSecondaryInt(GlobalIntegral* gq) override;
  //! \brief Returns the system quantity to be integrated by \a *this.
  GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const override;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
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

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] iGP Global integration point counter of first point in element
  bool finalizeElement(LocalIntegral& elmInt,
                       const FiniteElement& fe,
                       const TimeDomain& time, size_t iGP) override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time,
               const Vec3& X) const override;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Sets up the permeability matrix.
  //! \param[out] K \f$ n_{sd}\times n_{sd}\f$-matrix or its inverse
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] inverse If \e true, set up the inverse matrix instead
  bool formKmatrix(Matrix& K, const Vec3& X, bool inverse = false) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? 2*nsd+2 : 1; }

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

  //! \brief Returns order of time integration.
  int getOrder() const { return bdf.getActualOrder(); }

  //! \brief Set gravitational acceleration.
  void setGravity(double ga) { gacc = ga; }

  //! \brief Obtain integrand-type dependent solution norms
  virtual void getSolutionNorms(const SIMbase& sim, const Vector& solution,
                                double& dNorm,
                                double* dMax, size_t* iMax) const;

  //! \brief Returns pressure in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  virtual double pressure(const Vectors& eV,
                          const FiniteElement& fe,
                          size_t level) const;

  //! \brief Returns pressure gradient in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  virtual Vec3 pressureGradient(const Vectors& eV,
                                const FiniteElement& fe,
                                size_t level) const;

  //! \brief Helper for CoSTA.
  virtual void setParam(const std::string&, double) {}

  //! \brief Set material parameters.
  void setMaterial(DarcyMaterial& mat1) { mat = &mat1; }

  //! \brief Returns a const ref to material parameters.
  const DarcyMaterial& getMaterial() const { return *mat; }

  //! \brief Evaluates the darcy velocity in a point.
  bool evalDarcyVel (Vector& s, const Vectors& eV,
                     const FiniteElement& fe, const Vec3& X) const;

  //! \brief Returns gravitational acceleration.
  double getGravity() const { return gacc; }

  //! \brief Initializes and toggles the use of left-hand-side matrix buffers.
  void initLHSbuffers(size_t) override;

  //! \brief Enable/disable caching of element matrices.
  void lCache(bool enable) { useLCache = enable; }

  //! \brief Returns whether or not caching of element matrices is enabled.
  bool lCache() { return useLCache; }

protected:
  VecFunc*  bodyforce;    //!< Body force function
  VecFunc*  vflux;        //!< Flux function
  RealFunc* flux;         //!< Flux function
  RealFunc* source;       //!< Source function

  DarcyMaterial* mat = nullptr; //!< Material parameters

  GlobalIntegral* reacInt; //!< Reaction-forces-only integral
  TimeIntegration::BDF bdf; //!< BDF helper class

  double gacc = 9.81; //!< Gravity acceleration
  bool useLCache = true; //!< True to enable caching of element matrices
  bool reuseMats = false; //!< True to reuse matrices
  Matrices myKmats; //!< Cached element matrices

public:
  char extEner; //!< If \e true, external energy is to be computed
};


/*!
  \brief Class representing the integrand of Darcy energy norms.
*/

class DarcyNorm : public NormBase
{
public:
  //! \brief Enumeration of regular norm entries
  enum NormEntries {
    H1_Ph = 0,
    EXT_ENERGY,
    H1_Ch,
    H1_P,
    H1_E_Ph,
    H1_C,
    H1_E_Ch,
    TOTAL_NORM_E
  };

  //! \brief Enumeration of recovery norm entries
  enum RecoveryEntries {
    H1_Pr = 0,
    H1_Pr_Ph,
    H1_Cr,
    H1_Cr_Ch,
    TOTAL_NORM_REC,
    H1_E_Pr,
    H1_E_Cr,
    TOTAL_E_REC,
    EFF_REC_Ph,
    EFF_REC_Ch,
    EFF_REC_TOTAL
  };

  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical heat flux (optional)
  explicit DarcyNorm(Darcy& p, VecFunc* a = nullptr);
  //! \brief Empty destructor.
  virtual ~DarcyNorm();

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using NormBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite Element quantities
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt,  const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  bool finalizeElement(LocalIntegral& elmInt, const TimeDomain&,size_t) override;

  //! \brief Returns whether this norm has explicit boundary contributions.
  bool hasBoundaryTerms() const override { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

  //! \brief Projected quantities given as a field (recovery on a separate basis).
  //! \param field The field
  //! \param[in] idx The projection index
  void setProjectedFields(Fields* field, size_t idx) override;

  //! \brief Returns whether a norm quantity stores element contributions.
  bool hasElementContributions(size_t, size_t) const override;

protected:
  VecFunc* anasol; //!< Analytical heat flux
  std::vector<std::unique_ptr<Fields>> projFields; //!< Projected fields for recovery
};

#endif
