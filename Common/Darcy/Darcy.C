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
#include "Fields.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "SIMbase.h"
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
  bdf(torder)
{
  primsol.resize(1 + torder);

  tflux = nullptr;
  vflux = bodyforce = nullptr;
  flux = nullptr;
  reacInt = nullptr;
  extEner = false;
}


Darcy::~Darcy() = default;


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
  else if (tflux)
    return (*tflux)(X,normal)*normal;
  else
    return 0.0;
}


void Darcy::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  if (mode == SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.resize(1+bdf.getActualOrder());
}


void Darcy::setSecondaryInt (GlobalIntegral* gq)
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

  result->rhsOnly = neumann || this->reuseMats;
  result->withLHS = !neumann;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(nen);

  return result;
}


bool Darcy::initElement (const std::vector<int>& MNPC,
                         const FiniteElement& fe,
                         const Vec3& XC,
                         size_t nPt, LocalIntegral& elmInt)
{
  if (fe.iel > 0 && this->reuseMats)
  {
    size_t iel = fe.iel - 1;
    ElmMats* A = dynamic_cast<ElmMats*>(&elmInt);
    if (A && iel < this->myKmats.size() && !A->A.empty())
      A->A[0] = this->myKmats[iel];
  }

  return this->IntegrandBase::initElement(MNPC,fe,XC,nPt,elmInt);
}


bool Darcy::finalizeElement (LocalIntegral& elmInt,
                             const FiniteElement& fe,
                             const TimeDomain& time, size_t iGP)
{
  if (fe.iel > 0 && !this->reuseMats)
  {
    size_t iel = fe.iel - 1;
    ElmMats* A = dynamic_cast<ElmMats*>(&elmInt);
    if (A && iel < this->myKmats.size())
      this->myKmats[iel] = A->getNewtonMatrix();
  }

  return this->IntegrandBase::finalizeElement(elmInt,fe,time,iGP);
}


bool Darcy::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                     const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A.empty() && !this->reuseMats)
  {
    // Evaluate the hydraulic conductivity matrix at this point
    Matrix K;
    this->formKmatrix(K,X);
    WeakOps::LaplacianCoeff(elMat.A[pp], K, fe, 1.0/(this->getMaterial().rhow*gacc));
  }

  if (bodyforce)
  {
    // Integrate rhs contribution from body force
    Vec3 eperm = (*bodyforce)(X);
    Vec3 perm = this->getMaterial().getPermeability(X);
    for (size_t i = 0; i < nsd; i++)
      eperm[i] *= perm[i] / gacc;

    WeakOps::Divergence(elMat.b[pp], fe, eperm);
  }

  if (source)
    WeakOps::Source(elMat.b[pp], fe, (*source)(X));

  if (bdf.getActualOrder() > 0 && elmInt.vec.size() > 1 && pp == 0) {
    double p = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = this->pressure(elmInt.vec, fe, t);
      p -= val * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b[pp], fe, p);
    if (!elMat.A.empty() && !this->reuseMats)
      WeakOps::Mass(elMat.A[pp], fe, bdf[0] / time.dt);
  }

  if (m_mode == SIM::RHS_ONLY && !elmInt.vec.empty() && reacInt)
  {
    // Integrate the internal forces based on current solution
    Vector q;
    if (!this->evalDarcyVel(q,elmInt.vec,fe,X))
      return false;
    if (!fe.dNdX.multiply(q,elMat.b.front(),fe.detJxW,1.0)) // b += dNdX * q
      return false;
  }

  return true;
}


bool Darcy::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                     const Vec3& X, const Vec3& normal) const
{
  if (!flux && !vflux && !tflux)
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

  WeakOps::Source(elMat.b[pp], fe, qw / this->getMaterial().rhow);

  return true;
}


bool Darcy::formKmatrix (Matrix& K, const Vec3& X, bool inverse) const
{
  bool K_ok = true;
  Vec3 perm = this->getMaterial().getPermeability(X);

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
  if (!this->evalDarcyVel(s,eV,fe,X))
    return false;

  s.push_back(source ? (*source)(X) : 0.0);
  s.push_back(this->getMaterial().getPorosity(X));
  Vec3 perm = this->getMaterial().getPermeability(X);
  for (size_t i = 0; i < nsd; ++i)
    s.push_back(perm[i]);

  return true;
}


bool Darcy::evalDarcyVel (Vector& s, const Vectors& eV,
                          const FiniteElement& fe, const Vec3& X) const
{
  Vec3 dP = this->pressureGradient(eV, fe, 0);
  if (bodyforce)
    dP -= this->getMaterial().rhow * (*bodyforce)(X);

  Matrix K;
  this->formKmatrix(K,X);

  return K.multiply(Vector(dP.ptr(),nsd), s,
                    -1.0/(this->getMaterial().rhow*gacc));
}


std::string Darcy::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "pressure";

  return prefix + std::string(" pressure");
}


std::string Darcy::getField2Name (size_t i, const char* prefix) const
{
  if (i >= 2*nsd+2u) return "";

  if (nsd == 2 && i > 1)
   ++i;

  static const char* s[8] = {"v_x","v_y","v_z", "source", "porosity", "perm_x", "perm_y", "perm_z"};
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* Darcy::getNormIntegrand (AnaSol* asol) const
{
  return new DarcyNorm(*const_cast<Darcy*>(this),
                       asol ? asol->getScalarSecSol() : nullptr);
}


Vec3 Darcy::getBodyForce (const Vec3& X) const
{
  Vec3 result;
  if (bodyforce)
    result = (*bodyforce)(X);

  return result;
}


void Darcy::getSolutionNorms (const SIMbase& sim,
                              const Vector& solution,
                              double& dNorm,
                              double* dMax, size_t* iMax) const
{
  dNorm = sim.solutionNorms(solution, dMax, iMax);
}


Vec3 Darcy::pressureGradient (const Vectors& eV,
                              const FiniteElement& fe,
                              size_t level) const
{
  RealArray dP(nsd);
  fe.dNdX.multiply(eV.front(),dP,true);

  return dP;
}


double Darcy::pressure (const Vectors& eV,
                        const FiniteElement& fe,
                        size_t level) const
{
  return fe.N.dot(eV[level]);
}


void Darcy::initLHSbuffers (size_t nEl)
{
  if (!this->useLCache)
    return;

  if (nEl > 1) {
    this->myKmats.resize(nEl);
    this->reuseMats = false;
  } else if (nEl == 1)
    this->reuseMats = false;
  else if (nEl == 0 && !this->myKmats.empty())
    this->reuseMats = true;
}


DarcyNorm::DarcyNorm (Darcy& p, VecFunc* a) : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


DarcyNorm::~DarcyNorm() = default;


bool DarcyNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X) const
{
  Darcy& problem = static_cast<Darcy&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Kinv;
  if (!problem.formKmatrix(Kinv,X,true))
    return false; // Singular constitutive matrix

  double rgw = problem.getMaterial().rhow * problem.getGravity() * fe.detJxW;

  // Evaluate the finite element pressure field
  Vector dPh, dP, error;
  if (!problem.evalDarcyVel(dPh,pnorm.vec,fe,X))
    return false;

  // Integrate the energy norm a(p^h,p^h)
  pnorm[H1_Ph] += dPh.dot(Kinv*dPh)*rgw;
  // Evaluate the pressure field
  double p = problem.pressure(pnorm.vec, fe, 0);
  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[EXT_ENERGY] += problem.getPotential(X)*p*fe.detJxW;

  if (anasol)
  {
    // Evaluate the analytical velocity
    dP.fill((*anasol)(X).ptr(),fe.dNdX.cols());
    // Integrate the energy norm a(p,p)
    pnorm[H1_P] += dP.dot(Kinv*dP)*rgw;
    // Integrate the error in energy norm a(p-p^h,p-p^h)
    error = dP - dPh;
    double E = error.dot(Kinv*error)*rgw;
    pnorm[H1_E_Ph] += E;
    pnorm[TOTAL_NORM_E] += E;
  }

  size_t ip = this->getNoFields(1);
  size_t f = 2;
  for (const Vector& psol : pnorm.psol) {
    if (!projFields.empty() || !psol.empty())
    {
      Vector dPr(fe.dNdX.cols());
      if (!projFields.empty() && projFields[f-2]) {
        Vector vals;
        projFields[f-2]->valueFE(fe, vals);
        std::copy(vals.begin(), vals.begin()+fe.dNdX.cols(), dPr.begin());
      } else {
        // Evaluate projected pressure gradient
        for (size_t j = 0; j < fe.dNdX.cols(); j++)
          dPr[j] = psol.dot(fe.N,j,nrcmp);
      }

      // Integrate the energy norm a(p^r,p^r)
      pnorm[ip+H1_Pr] += dPr.dot(Kinv*dPr)*rgw;

      // Integrate the estimated error in energy norm a(p^r-p^h,p^r-p^h)
      error = dPr - dPh;
      double E = error.dot(Kinv*error)*rgw;
      pnorm[ip+H1_Pr_Ph] += E;
      pnorm[ip+TOTAL_NORM_REC] += E;

      if (anasol)
      {
        // Integrate the error in the projected solution a(p-p^r,p-p^r)
        error = dP - dPr;
        E = error.dot(Kinv*error)*rgw;
        pnorm[ip+H1_E_Pr] += E;
        pnorm[ip+TOTAL_E_REC] += E;
      }
    }
    ip += this->getNoFields(f++);
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

  size_t ip = this->getNoFields(1);
  for (size_t k = 0; k < pnorm.psol.size() && ip < pnorm.size(); ++k) {
    pnorm[ip+EFF_REC_Ph] = sqrt(pnorm[ip+H1_Pr_Ph] / pnorm[H1_E_Ph]);
    pnorm[ip+EFF_REC_Ch] = sqrt(pnorm[ip+H1_Cr_Ch] / pnorm[H1_E_Ch]);
    pnorm[ip+EFF_REC_TOTAL] = sqrt((pnorm[ip+H1_Pr_Ph] + pnorm[ip+H1_Cr_Ch]) /
                                   (pnorm[H1_E_Ph] + pnorm[H1_E_Ch]));
    ip += this->getNoFields(k+2);
  }

  return true;
}


size_t DarcyNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return group == 1 ? 8 : 11;
}


std::string DarcyNorm::getName (size_t i, size_t j, const char* prefix) const
{
  static const char* s[8] = {
    "a(p^h,p^h)^0.5",
    "(h,p^h)^0.5",
    "a(c^h,c^h)^0.5",
    "a(p,p)^0.5",
    "a(e,e)^0.5, e=p-p^h",
    "a(c,c)^0.5",
    "a(e,e)^0.5, e=c-c^h",
    "a(e,e)^0.5, e=(p,c)-(p,c)^h"
  };

  static const char* r[11] = {
    "a(p^r,p^r)^0.5",
    "a(e,e)^0.5, e=p^r-p^h",
    "a(c^r,c^r)^0.5",
    "a(e,e)^0.5, e=c^r-c^h",
    "|e|, (p,c)^r-(p,c)^h",
    "a(e,e)^0.5, e=p-p^r",
    "a(e,e)^0.5, e=c-c^r",
    "|e|, e=(p,c)^r-(p,c)",
    "eta^p",
    "eta^c",
    "eta^tot"
  };

  const char** n = i > 1 ? r : s;

  if (!prefix)
    return n[j-1];

  return prefix + std::string(" ") + n[j-1];
}


void DarcyNorm::setProjectedFields (Fields* field, size_t idx)
{
  if (idx >= projFields.size())
    projFields.resize(idx+1);
  projFields[idx].reset(field);
}


bool DarcyNorm::hasElementContributions(size_t i, size_t j) const
{
  if (i == 1) {
    if (j == 2)
      return false;
    if (!anasol && j > 3)
     return false;
   } else {
    if (!anasol && j > 5)
      return false;
  }

  return true;
}
