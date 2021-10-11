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

#include "Function.h"
#include "FunctionSum.h"

#include <vector>


class Vec3;


/*!
 \brief Primary solution for Darcy problem on a L-shape domain.
*/

class LshapeDarcy : public RealFunc
{
public:
  //! \brief Empty constructor.
  LshapeDarcy() {}
  //! \brief Empty destructor.
  virtual ~LshapeDarcy() {}

protected:
  //! \brief Evaluates the analytic pressure field at the point \a X.
  Real evaluate(const Vec3& X) const override;
};


/*!
 \brief Secondary solution for Darcy problem on a L-shape domain.
*/

class LshapeDarcyVelocity : public VecFunc
{
public:
  //! \brief Empty constructor.
  LshapeDarcyVelocity() {}
  //! \brief Empty destructor.
  virtual ~LshapeDarcyVelocity() {}

protected:
  //! \brief Evaluates the analytic vecocity vector at the point \a X.
  Vec3 evaluate(const Vec3& X) const override;
};


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

  //! \brief Parse functions from a string.
  //! \param input The string to parse
  bool parse(const char* input);

  //! \brief Empty constructor.
  virtual ~DiracSum() {}

protected:
  //! \brief Evaluates the function in a point.
  //! \param X Coordinates of point to evaluate in
  Real evaluate(const Vec3& X) const override
  { return this->FunctionSum::getValue(X).front(); }

  double pointTol; //!< Interval around point associated with functions
  int myDim; //!< Dimensionality of world
};

#endif
