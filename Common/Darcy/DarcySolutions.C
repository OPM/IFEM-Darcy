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

#include "ASMbase.h"
#include "ExprFunctions.h"
#include "IFEM.h"
#include "LogStream.h"
#include "SIMbase.h"
#include "StringUtils.h"
#include "Vec3.h"

#include <array>
#include <cmath>
#include <sstream>
#include <string>


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
    Real x, y = 0.0, z = 0.0;
    s2 >> x;
    if (myDim > 1)
      s2 >> y;
    if (myDim > 2)
      s2 >> z;
    std::string value;
    s2 >> value;
    IFEM::cout << "\n\t\tDirac(" << x;
    if (myDim > 1)
      IFEM::cout << ", " << y;
    if (myDim > 2)
      IFEM::cout << ", " << z;
    IFEM::cout << ") = " << value;

    // Computed with quadrature - needed to normalize integral to 1
    double normfactor =
        (myDim == 1) ? 0.44399381616807944 :
        (myDim == 2) ? 0.46651239317833007 :
                       0.44108888727660440;
    std::stringstream s3;
    s3.precision(17);
    s3 << "r2=pow(x-" << x << ",2)";
    if (myDim > 1)
      s3 << "+pow(y-" << y << ",2)";
    if (myDim > 2)
      s3 << "+pow(z-" << z << ",2)";
    s3 << "; r=sqrt(r2); "
       << "eps=" << pointTol << "; "
       << "n=" << normfactor << "; "
       << "eps2=eps*eps; "
       << "epsd=pow(eps," << myDim << "); "
       << "(" << value << ")"
       << "*below(r,eps)*(1/(epsd*n))*exp(-eps2/(eps2-r2))";
    m_funcs.push_back(std::make_unique<EvalFunction>(s3.str().c_str()));
    this->add(m_funcs.back().get());
    ok = true;
  }

  IFEM::cout << std::endl;
  return ok;
}


void DiracSum::setParam (const std::string& name, double value)
{
  for (const auto& func : m_funcs)
    static_cast<EvalFunction*>(func.get())->setParam(name,value);
}


bool ElementSum::parse (const char* input, const SIMbase& sim)
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
    std::array<double,3> p{0.0,0.0,0.0};
    s2 >> p[0];
    if (myDim > 1)
      s2 >> p[1];
    if (myDim > 2)
      s2 >> p[2];
    std::string value;
    s2 >> value;
    int patch = 1;
    s2 >> patch;

    IFEM::cout << "\n\t\tElement(" << p[0];
    if (myDim > 1)
      IFEM::cout << ", " << p[1];
    if (myDim > 2)
      IFEM::cout << ", " << p[2];
    IFEM::cout << ", " << patch << ") = " << value;

    ASMbase* pch = sim.getPatch(patch);
    if (!patch) {
      std::cerr << "** ElementSum: No patch " << patch << ".\n";
      continue;
    }

    int iel = pch->findElementContaining(p.data());
    if (iel <= 0) {
      std::cerr << "** ElementSum: Failed to locate element\n";
      continue;
    }

    Matrix X;
    if (!pch->getElementCoordinates(X, iel)) {
      std::cerr << "** ElementSum: Failed to obtain element coordinates\n";
      continue;
    }

    std::stringstream s3;
    s3 << "above(x+0.001," << X(1,1) << ")*"
       << "below(x-0.001," << X(1,2) << ")";
    if (myDim > 1)
       s3 << "*above(y+0.001," << X(2,1) << ")*"
          << "below(y-0.001," << X(2,3) << ")";
    if (myDim > 2)
       s3 << "*above(z+0.001," << X(3,1) << ")*"
          << "below(z-0.001," << X(3,8) << ")";
    s3 << "*" << value;
    IFEM::cout << " -> " << s3.str();
    m_funcs.push_back(std::make_unique<EvalFunction>(s3.str().c_str()));
    this->add(m_funcs.back().get());
    ok = true;
  }

  IFEM::cout << std::endl;
  return ok;
}


void ElementSum::setParam (const std::string& name, double value)
{
  for (const auto& func : m_funcs)
    static_cast<EvalFunction*>(func.get())->setParam(name,value);
}
