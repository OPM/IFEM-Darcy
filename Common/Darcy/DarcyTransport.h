// $Id$
//==============================================================================
//!
//! \file DarcyTransport.h
//!
//! \date Oct 20 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy flow with tracer transport.
//!
//==============================================================================

#ifndef _DARCY_TRANSPORT_H_
#define _DARCY_TRANSPORT_H_

#include "Darcy.h"

#include <memory>


/*!
  \brief Class representing the integrand of the Darcy transport problem.
*/

class DarcyTransport : public Darcy
{
public:
  //! \brief The constructor initializes all pointers to zero.
  explicit DarcyTransport(unsigned short int n, int torder = 0);
  //! \brief Empty destructor.
  ~DarcyTransport() override;

  //! \brief Defines the concentration source function.
  void setCSource(std::unique_ptr<RealFunc> s) override { sourceC = std::move(s); }

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt) override;

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

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  bool evalSol (Vector& s,
                const FiniteElement& fe,
                const Vec3& X,
                const std::vector<int>& MNPC) const override;

  using IntegrandBase::finalizeElement;

  //! \brief Finalizes the element quantities after the numerical integration.
  bool finalizeElement(LocalIntegral& A) override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? 3*nsd+3 : 2; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Index for field
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

  //! \brief Returns pressure in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  double pressure(const Vectors& eV,
                  const FiniteElement& fe,
                  size_t level) const override;

  //! \brief Returns concentration gradient in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  Vec3 concentrationGradient(const Vectors& eV,
                             const FiniteElement& fe,
                             size_t level) const;

  //! \brief Returns pressure gradient in a point.
  //! \param eV Element vectors
  //! \param fe Finite element data at current point
  //! \param level Time level to evaluate at
  Vec3 pressureGradient(const Vectors& eV,
                        const FiniteElement& fe,
                        size_t level) const override;

protected:
  std::unique_ptr<RealFunc> sourceC; //!< Concentration source function
};


/*!
  \brief Class representing the integrand of Darcy transport energy norms.
*/

class DarcyTransportNorm : public DarcyNorm
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical pressure gradient (optional)
  //! \param[in] c The analytical concentration gradient (optional)
  explicit DarcyTransportNorm(DarcyTransport& p,
                              VecFunc* a = nullptr,
                              VecFunc* c = nullptr);

  //! \brief Empty destructor.
  virtual ~DarcyTransportNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

private:
  VecFunc* anac; //!< Analytical concentration gradient
};

#endif
