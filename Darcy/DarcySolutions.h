// $Id$
//==============================================================================
//!
//! \file DarcySolutions.h
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Analytic solutions for Darcy problems.
//!
//==============================================================================

#ifndef _DARCY_SOLUTIONS_H
#define _DARCY_SOLUTIONS_H

#include "FunctionSum.h"

#include <memory>
#include <vector>

class SIMbase;


/*!
 \brief A sum of dirac functions as a source.
*/

class DiracSum : public FunctionSum, public RealFunc
{
public:
  //! \brief Constructor.
  //! \param tol Interval around point associated with functions
  //! \param dim Dimensionality of world
  DiracSum(double tol, int dim) : pointTol(tol), myDim(dim) {}

  //! \brief Empty constructor.
  virtual ~DiracSum() {}

  //! \brief Parse functions from a string.
  //! \param input The string to parse
  bool parse(const char* input);

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value) override;

protected:
  //! \brief Evaluates the function in a point.
  //! \param X Coordinates of point to evaluate in
  Real evaluate(const Vec3& X) const override
  { return this->FunctionSum::getValue(X).front(); }

  double pointTol; //!< Interval around point associated with functions
  int myDim; //!< Dimensionality of world

  std::vector<std::unique_ptr<RealFunc>> m_funcs; //!< Vector of pointers to functions
};


/*!
 \brief A sum of single-element sources.
*/

class ElementSum : public FunctionSum, public RealFunc
{
public:
  //! \brief Constructor.
  //! \param dim Dimensionality of world
  ElementSum(int dim) : myDim(dim) {}

  //! \brief Empty constructor.
  virtual ~ElementSum() {}

  //! \brief Parse functions from a string.
  //! \param input The string to parse
  //! \param sim Simulator with elements information
  bool parse(const char* input, const SIMbase& sim);

  //! \brief Set an additional parameter in the function.
  void setParam(const std::string& name, double value) override;

protected:
  //! \brief Evaluates the function in a point.
  //! \param X Coordinates of point to evaluate in
  Real evaluate(const Vec3& X) const override
  { return this->FunctionSum::getValue(X).front(); }

  int myDim; //!< Dimensionality of world
  std::vector<std::unique_ptr<RealFunc>> m_funcs; //!< Vector of pointers to functions
};

#endif
