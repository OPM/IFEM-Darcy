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
    virtual double evaluate(const Vec3& X) const;
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
    virtual Vec3 evaluate(const Vec3& X) const;
};

/*!
 \brief Primary solution for a wavefront Darcy problem.
*/

class Wavefront : public RealFunc
{
  public:
    //! \brief Empty constructor
    Wavefront() {}
    //! \brief Empty destructor
    virtual ~Wavefront() {}
  protected:
    //! \brief Evaluates the analytical pressure field at the point \a X.
    virtual double evaluate(const Vec3& X) const;
};


/*!
 \brief Secondary solution for a wavefront Darcy problem.
*/

class WavefrontVelocity : public VecFunc
{
  public:
    //! \brief Empty constructor
    WavefrontVelocity() {}
    //! \brief Empty destructor
    virtual ~WavefrontVelocity() {}
  protected:
    //! \brief Evaluates the analytical pressure field at the point \a X.
    virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
 \brief Source function for a wavefront Darcy problem.
*/

class WavefrontSource : public RealFunc
{
  public:
    //! \brief Empty constructor
    WavefrontSource() {}
    //! \brief Empty destructor
    virtual ~WavefrontSource() {}
  protected:
    //! \brief Evaluates the analytical source at the point \a X.
    virtual double evaluate(const Vec3& X) const;
};

#endif
