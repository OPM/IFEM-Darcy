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

#ifndef DARCY_H_
#define DARCY_H_

#include "IntegrandBase.h"
#include "ElmMats.h"
#include "EqualOrderOperators.h"
#include "Vec3.h"

class DarcyNorm;
class RealFunc;
class VecFunc;


/*!
 *  \brief Class representing the integrand of the Darcy problem.
 */
class Darcy : public IntegrandBase
{
public:
  using WeakOps = EqualOrderOperators::Weak; //!< Convenience rename
  //! \brief Default constructor.
  explicit Darcy(unsigned short int n);
  //! \brief Empty destructor.
  virtual ~Darcy() {}

  //! \brief Set permeability function.
  void setPermValues(VecFunc* permv) { permvalues = permv; }

  //! \brief Returns permeability in a point
  Vec3 getPermeability(const Vec3& X) const;

  //! \brief Set permeability field.
  void setPermField(RealFunc* permf) { permeability = permf; }

  //! \brief Set body force vector.
  void setBodyForce(VecFunc* body) { bodyforce = body; }

  //! \brief Returns body forces in a point
  Vec3 getBodyForce(const Vec3& X) const;

  //! \brief Defines the source function.
  void setSource(RealFunc* src) { source = src; }

  //! \brief Defines a scalar flux function.
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Defines a vectorial flux function.
  void setFlux(VecFunc* vf) { vflux = vf; }

  //! \brief Evaluates the boundary fluid flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& normal) const;
  //! \brief Evaluates the potential source (if any) at specified point.
  double getPotential(const Vec3& X) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Sets up the permeability matrix
  //! \param[out] K \f$ nsd\times nsd\f$-matrix or its inverse
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] inverse If \e true, set up the inverse matrix instead
  bool formKmatrix(Matrix& K, const Vec3& X, bool inverse = false) const;

  using IntegrandBase::evalSol;
  //! \brief Evaluated the secondary solution at a result point
  //! param[out] s Array of solution field values at current point
  //! param[in] fe Finite element data at current point
  //! param[in] X Cartesian coordinates of current point
  //! param(in) MNPC Nodal point correspondence for the basis function values
  bool evalSol(Vector& s, const FiniteElement& fe,
               const Vec3& X, const std::vector<int>& MNPC) const override;

  //! \brief Evaluates the finite element (FE) solution at an integration point
  //! \param[out] s The FE solution values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol(Vector& s, const Vector& eV,
               const Matrix& dNdX, const Vec3& X) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld = 2) const override { return fld > 1 ? nsd : 1; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix = 0) const override;

  //! \brief Returns the name of a secondary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix = 0) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol = 0) const override;

private:
  VecFunc* permvalues;    //!< Permeability function.
  RealFunc* permeability; //!< Permeability field function.

protected:
  RealFunc* flux;         //!< Flux function.
  RealFunc* source;       //!< Source function.
  VecFunc* vflux;         //!< Flux function.

public:
  const double rhow;      //!< Density of fluid.
  const double gacc;      //!< Gravity acceleration.
  VecFunc* bodyforce;     //!< Body force function.
};


/*!
 *   \brief Class representing the integrand of Darcy energy norms.
 */

class DarcyNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical heat flux (optional)
  DarcyNorm(Darcy& p, VecFunc* a = 0);
  //! \brief Empty destructor.
  virtual ~DarcyNorm() {}

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
  size_t getNoFields(int group = 0) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

  //! \brief Returns whether a norm quantity stores element contributions.
  bool hasElementContributions(size_t i, size_t j) const override;

private:
  VecFunc* anasol; //!< Analytical heat flux

};

#endif /* DARCY_H_ */
