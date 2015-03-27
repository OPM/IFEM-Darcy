// $Id$
//==============================================================================
//!
//! \file DarcySolutions.C
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Analytic solutions for Darcy problems.
//!
//==============================================================================

#include "DarcySolutions.h"
#include "Vec3.h"

double LshapeDarcy::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  double r2 = x*x + y*y;
  double theta = atan2(y,x);
  if (theta <= 0) theta += 2*pi;

  return pow(r2,1.0/3)*(sin((2.0*theta-pi)/3));
}

Vec3 LshapeDarcyVelocity::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  double r2 = x*x + y*y;
  double theta = atan2(y,x);
  if (theta <= 0) theta += 2*pi;
  if (r2 < 1e-16) {  // truncate the singularity to avoid NaN values
    r2 = 1e-16;
    theta = 2*pi;
  }

  Vec3 velocity;

  velocity.x = (2.0/3) * (cos(2.0/3*theta + pi/6)*x + sin(2.0/3*theta + pi/6)*y) / pow(r2, 2.0/3);
  velocity.y = (-2.0/3) * (sin(2.0/3*theta - pi/3)*y + sin(2.0/3*theta + pi/6)*x) / pow(r2, 2.0/3);

  return velocity;
}

double Wavefront::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;

  double a = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));

  return atan(50.0*(a-0.25));
}

Vec3 WavefrontVelocity::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;

  double a = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
  double b = 1.0 + 2500.0*(a-0.25)*(a-0.25);
  double c = 10.0; // Degree of anisotropy

  Vec3 velocity;
  velocity.x = -50.0*c*(x-0.5)/(a*b);
  velocity.y = -50.0*(y-0.5)/(a*b);

  return velocity;
}

double WavefrontSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;

  double a = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
  double b = 1.0 + 2500.0*(a-0.25)*(a-0.25);
  double c = 10.0; // Degree of anisotropy

  double f1 = 250000.0*(a-0.25)/(a*a*b*b);
  double f2 = 50.0/(a*a*a*b);
  double f3 = 50.0*(c+1)/(a*b);

  double f = (f1+f2)*(c*(x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) - f3;

  return f;
}
