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


/*!
  \brief Class representing the integrand of the transport correction problem.
*/

class DarcyTransportCorr : public IntegrandBase
{
public:
  //! \brief Enum for element level right-hand-side vectors
  enum ResidualVectors
  {
    Fq = 1, Fl = 2, Fg = 3, NVEC = 4
  };

  //! \brief Enum for element level left-hand-side matrices
  enum TangentMatrices
  {
    qq   = 1, ql = 4,
    lq   = 6,         lg = 7,
    NMAT = 10
  };

  //! \brief Default constructor.
  explicit DarcyTransportCorr(unsigned short int n = 3);
  //! \brief Default destructor.
  ~DarcyTransportCorr() override;

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

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element quantities after the numerical integration.
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents.
  bool finalizeElement(LocalIntegral& elmInt) override;

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  bool evalSol(Vector& s,
               const FiniteElement& fe, const Vec3& X,
               const std::vector<int>& MNPC) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? nsd+3 : nsd; }

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

  //! \brief Set observed concentration.
  void setObservedConcentration(std::unique_ptr<RealFunc> obs_c);

  //! \brief Set input source function.
  void setInputSource(std::unique_ptr<RealFunc> f);

  //! \brief Set input Darcy velocity.
  void setInputVelocity(std::unique_ptr<VecFunc> inp_q);

  //! \brief Set penalty parameter.
  void setMassPenaltyParam(const double nmu) { alpha = nmu;}

  //! \brief Set penalty parameter.
  void setTransportPenaltyParam(const double nmu) { beta = nmu;}

private:
  double residual(const Vec3& X) const;

  double alpha = 1.0e6;  //!< Mass penalty parameter
  double beta  = 1.0e6;  //!< Transport penalty parameter
  std::unique_ptr<VecFunc>  input_q;      //!< Input Darcy velocity
  std::unique_ptr<RealFunc> input_source; //!< Input source
  std::unique_ptr<RealFunc> observed_C;   //!< Observed tracer concentration
};

#endif
