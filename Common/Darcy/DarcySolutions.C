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

#include "ExprFunctions.h"
#include "IFEM.h"
#include "StringUtils.h"
#include "Vec3.h"

#include <sstream>

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


bool DiracSum::parse (const char* input)
{
  std::stringstream str;
  std::string val = input;
  replaceAll(val, "\\", "\n");
  str << val;
  bool ok = false;
  while (str.good()) {
    char temp[1024];
    str.getline(temp, 1024);
    if (temp[0] == '#' || temp[0] == 0)
      continue;
    std::stringstream s2(temp);
    Real x, y = 0.0, z = 0.0, value;
    s2 >> x;
    if (myDim > 1)
      s2 >> y;
    if (myDim > 2)
      s2 >> z;
    s2 >> value;
    IFEM::cout << "\n\t\tDirac(" << x;
    if (myDim > 1)
      IFEM::cout << ", " << y;
    if (myDim > 2)
      IFEM::cout << ", " << z;
    IFEM::cout << ") = " << value;
    std::stringstream s3;
    s3 << "r2=pow(x-" << x << ",2)";
    if (myDim > 1)
      s3 << "+pow(y-" << y << ",2)";
    if (myDim > 2)
      s3 << "+pow(z-" << z << ",2)";
     s3 << "; r=sqrt(r2); if(below(r,"
       << pointTol << ")," << value << ",0.0)";
    EvalFunction* e = new EvalFunction(s3.str().c_str());
    this->add(e);
    ok = true;
  }

  IFEM::cout << std::endl;
  return ok;
}
