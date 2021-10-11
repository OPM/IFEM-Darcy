// $Id$
//==============================================================================
//!
//! \file Darcy.C
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for Darcy flow problems.
//!
//==============================================================================

#include "Darcy.h"

#include "AnaSol.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <cmath>
#include <ext/alloc_traits.h>
#include <iostream>
#include <memory>
#include <vector>


Darcy::Darcy (unsigned short int n, int torder) :
  IntegrandBase(n),
  bdf(torder),
  rhow(1.0),
  gacc(9.81)
{
  primsol.resize(1 + torder);

  permvalues = vflux = bodyforce = nullptr;
  permeability = flux = source = nullptr;
  reacInt = nullptr;
  extEner = false;
}


double Darcy::getPotential (const Vec3& X) const
{
  return source ? (*source)(X) : 0.0;
}


double Darcy::getFlux (const Vec3& X, const Vec3& normal) const
{
  if (flux)
    return (*flux)(X);
  else if (vflux)
    return (*vflux)(X)*normal;
  else
    return 0.0;
}


void Darcy::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;

  primsol.resize(1 + bdf.getActualOrder());
}


void Darcy::setReactionIntegral (GlobalIntegral* gq)
{
  delete reacInt;
  reacInt = gq;
}


GlobalIntegral& Darcy::getGlobalInt (GlobalIntegral* gq) const
{
  if (m_mode == SIM::RHS_ONLY && reacInt)
    return *reacInt;

  return this->IntegrandBase::getGlobalInt(gq);
}


LocalIntegral* Darcy::getLocalIntegral (size_t nen, size_t, bool neumann) const
{
  ElmMats* result = new ElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(nen);

  return result;
}


bool Darcy::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                     const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A.empty())
  {
    // Evaluate the hydraulic conductivity matrix at this point
    Matrix K;
    this->formKmatrix(K,X);
    WeakOps::LaplacianCoeff(elMat.A.front(), K, fe, 1.0/(rhow*gacc));
  }

  if (bodyforce)
  {
    // Integrate rhs contribution from body force
    Vec3 eperm = (*bodyforce)(X);
    if (permvalues)
    {
      // Read components of permeability matrix
      Vec3 perm = (*permvalues)(X);
      for (size_t i = 0; i < nsd; i++)
        eperm[i] *= perm[i]/gacc;
    }
    else if (permeability)
      // Read point-wise permeabilities from a random field
      eperm *= (*permeability)(X)/gacc;

    WeakOps::Divergence(elMat.b.front(), fe, eperm);
  }

  if (source)
    WeakOps::Source(elMat.b.front(), fe, (*source)(X));

  if (bdf.getActualOrder() > 0 && elmInt.vec.size() > 1) {
    double p = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = fe.N.dot(elmInt.vec[t]);
      p -= val * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b.front(), fe, p);
    WeakOps::Mass(elMat.A.front(), fe, bdf[0] / time.dt);
  }

  if (m_mode == SIM::RHS_ONLY && !elmInt.vec.empty() && reacInt)
  {
    // Integrate the internal forces based on current solution
    Vector q;
    if (!this->evalSol2(q,elmInt.vec,fe,X))
      return false;
    if (!fe.dNdX.multiply(q,elMat.b.front(),fe.detJxW,1.0)) // b += dNdX * q
      return false;
  }

  return true;
}


bool Darcy::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                     const Vec3& X, const Vec3& normal) const
{
  if (!flux && !vflux)
  {
    std::cerr <<" *** Darcy::evalBou: No fluxes."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty())
  {
    std::cerr <<" *** Darcy::evalBou: No load vector."<< std::endl;
    return false;
  }

  double qw = -this->getFlux(X,normal);

  WeakOps::Source(elMat.b.front(), fe, qw/rhow);

  return true;
}


bool Darcy::formKmatrix (Matrix& K, const Vec3& X, bool inverse) const
{
  bool K_ok = true;
  Vec3 perm = this->getPermeability(X);

  K.resize(nsd,nsd,true);
  for (int i = 1; i <= nsd; i++)
    if (!inverse)
      K(i,i) = perm[i-1];
    else if (perm[i-1] != 0.0)
      K(i,i) = 1.0 / perm[i-1];
    else
      K_ok = false;

  return K_ok;
}


bool Darcy::evalSol2 (Vector& s, const Vectors& eV,
                      const FiniteElement& fe, const Vec3& X) const
{
  // Evaluate the pressure gradient
  RealArray temp;
  if (eV.empty() || !fe.dNdX.multiply(eV.front(),temp,true))
  {
    std::cerr <<" *** Darcy::evalSol: Invalid solution vector.\n"
              <<"     size(eV) = "<< (eV.empty() ? 0 : eV.front().size())
              <<" size(dNdX) = "<< fe.dNdX.rows() <<","<< fe.dNdX.cols()
              << std::endl;
    return false;
  }

  if (bodyforce)
  {
    Vec3 b = (*bodyforce)(X);
    for (size_t i = 0; i < nsd; i++)
      temp[i] -= rhow*b[i];
  }

  Matrix K;
  this->formKmatrix(K,X);

  return K.multiply(temp,s,-1.0/(rhow*gacc));
}


std::string Darcy::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "pressure";

  return prefix + std::string(" pressure");
}


std::string Darcy::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return "";

  static const char* s[3] = {"v_x","v_y","v_z"};
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* Darcy::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new DarcyNorm(*const_cast<Darcy*>(this),asol->getScalarSecSol());
  else
    return new DarcyNorm(*const_cast<Darcy*>(this));
}


Vec3 Darcy::getPermeability (const Vec3& X) const
{
  Vec3 result;
  if (permvalues)
    result = (*permvalues)(X);
  else if (permeability)
    result = (*permeability)(X);

  return result;
}


Vec3 Darcy::getBodyForce (const Vec3& X) const
{
  Vec3 result;
  if (bodyforce)
    result = (*bodyforce)(X);

  return result;
}


DarcyNorm::DarcyNorm (Darcy& p, VecFunc* a) : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


bool DarcyNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X) const
{
  Darcy& problem = static_cast<Darcy&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Kinv;
  if (!problem.formKmatrix(Kinv,X,true))
    return false; // Singular constitutive matrix

  double rgw = problem.rhow * problem.gacc * fe.detJxW;

  // Evaluate the finite element pressure field
  Vector sigmah, sigma, error;
  if (!problem.evalSol2(sigmah,pnorm.vec,fe,X))
    return false;

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += sigmah.dot(Kinv*sigmah)*rgw;
  // Evaluate the pressure field
  double p = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[1] += problem.getPotential(X)*p*fe.detJxW;

  size_t ip = 2;
  if (anasol)
  {
    // Evaluate the analytical velocity
    sigma.fill((*anasol)(X).ptr(),nrcmp);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Kinv*sigma)*rgw;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Kinv*error)*rgw;
  }

  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate projected heat flux field
      Vector sigmar(nrcmp);
      for (size_t j = 0; j < nrcmp; j++)
        sigmar[j] = psol.dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Kinv*sigmar)*rgw;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Kinv*error)*rgw;

      if (anasol)
      {
        // Integrate the error in the projected solution a(u-u^r,u-u^r)
        error = sigma - sigmar;
        pnorm[ip++] += error.dot(Kinv*error)*rgw;
        ip++; // Make room for the local effectivity index here
      }
    }

  return true;
}


bool DarcyNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X, const Vec3& normal) const
{
  const Darcy& problem = static_cast<const Darcy&>(myProblem);
  if (!problem.extEner) return true;

  // Evaluate the surface heat flux
  double h = problem.getFlux(X,normal);
  // Evaluate the temperature field
  double u = elmInt.vec.front().dot(fe.N);

  // Integrate the external energy (h,u^h)
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  pnorm[1] += h*u*fe.detJxW;
  return true;
}


bool DarcyNorm::finalizeElement (LocalIntegral& elmInt,
                                 const TimeDomain&, size_t)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 7; ip < pnorm.size(); ip += 4)
    pnorm[ip] = sqrt(pnorm[ip-2] / pnorm[3]);

  return true;
}


size_t DarcyNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return anasol ? 4 : 2;
}


std::string DarcyNorm::getName (size_t i, size_t j, const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[8] = {
    "a(u^h,u^h)^0.5",
    "(h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index"
  };

  size_t k = i > 1 ? j+3 : j-1;

  if (!prefix)
    return s[k];

  return prefix + std::string(" ") + s[k];
}


bool DarcyNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}
