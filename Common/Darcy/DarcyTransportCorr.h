// $Id$
//==============================================================================
//!
//! \file DarcyTransportCorr.h
//!
//! \date Jan 14 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy flow with tracer transport.
//!
//==============================================================================

#ifndef _DARCY_TRANSPORT_CORR_H_
#define _DARCY_TRANSPORT_CORR_H_

#include "IntegrandBase.h"

#include <memory>

class RealFunc;
class Tensor;


/*!
  \brief Class representing the integrand of the transport correction problem.
*/

class DarcyTransportCorr : public IntegrandBase
{
public:
  //! \brief Enum for element level right-hand-side vectors
  enum ResidualVectors
  {
    Fq = 1, Fl = 2, Fm = 3, NVEC = 4
  };

  //! \brief Enum for element level left-hand-side matrices.
  //!
  //! \details The matrix layout for 3x3 blocks is as follows:
  //! \code
  //!     1 4 5
  //!     6 2 8
  //!     7 9 3
  //! \endcode

  enum TangentMatrices
  {
    qq = 1, ql = 4, qm = 5, NMAT = 10
  };

  //! \brief Default constructor.
  explicit DarcyTransportCorr(unsigned short int n = 3);
  //! \brief Default destructor.
  ~DarcyTransportCorr() override;

  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  LocalIntegral* getLocalIntegral(size_t nen,
                                  size_t, bool) const override;

  //! \brief Returns a local integral container for the given element
  //! \param[in] nen Number of nodes on element
  LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                  size_t, bool) const override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  using IntegrandBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  bool finalizeElement(LocalIntegral& elmInt) override;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? 5+nsd : nsd; }

  //! \brief Filters a result components for output.
  bool suppressOutput(size_t i, ASM::ResultClass type) const override;

  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Field index
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

  //! \brief Evaluates the FE solution.
  Vec3 evalSol(const Vector& eV, const Vector& N) const;
  //! \brief Evaluates the FE solution gradient.
  Tensor evalGrd(const Vector& eV, const Matrix& dNdX) const;

  //! \brief Evaluates the input velocity function.
  Vec3 evalInput(const Vec3& X) const;
  //! \brief Evaluates the divergence of the input velocity function.
  double divInput(const Vec3& X) const;
  //! \brief Evaluates the tracer concentration solution.
  double evalTracer(const Vec3& X) const;

  //! \brief Evaluates the residual using the input velocity.
  double residual(const Vec3& X) const;

  //! \brief Evaluates the residual in the FE solution.
  double residual(const Vector& eV,
                  const Vector& N, const Matrix& dNdX, const Vec3& X) const;

private:
  double alpha = 1.0e6;  //!< Mass penalty parameter
  double beta  = 1.0e6;  //!< Transport penalty parameter
  double eps   = 1.0e-6; //!< Division by zero tolerance in mass-term scaling
  std::unique_ptr<VecFunc>  input_q;      //!< Input Darcy velocity
  std::unique_ptr<RealFunc> input_source; //!< Input source
  std::unique_ptr<RealFunc> observed_C;   //!< Observed tracer concentration
};


/*!
  \brief Class representing the integrand of Darcy transport correction norms.
*/

class DarcyTCNorm : public NormBase
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit DarcyTCNorm(DarcyTransportCorr& p): NormBase(p) {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

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
