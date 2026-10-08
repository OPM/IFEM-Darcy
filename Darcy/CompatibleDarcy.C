// $Id$
//==============================================================================
//!
//! \file CompatibleDarcy.C
//!
//! \date Sep 21 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand for Darcy flow with a div-compatible basis.
//!
//==============================================================================

#include "CompatibleDarcy.h"
#include "DarcyMaterial.h"

#include "AnaSol.h"
#include "BlockElmMats.h"
#include "CompatibleOperators.h"
#include "ElmNorm.h"
#include "EqualOrderOperators.h"
#include "Fields.h"
#include "FiniteElement.h"
#include "Function.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3Oper.h"


CompatibleDarcy::CompatibleDarcy (unsigned short int n, int torder)
  : Darcy(n,torder)
{
  npv = n+1;
  Darcy::pp = pp;
  scalarBasis = n+1;
}


LocalIntegral* CompatibleDarcy::getLocalIntegral (const UintVec& nen,
                                                  size_t, bool neumann) const
{
  const size_t nBlock = 4;
  const size_t nBasis = 4;
  BlockElmMats* result = new BlockElmMats(nBlock,nBasis);

  size_t i = 0;
  result->resize(NMAT, NVEC);
  result->redim(qxqx, nen[i++], 1, 1);
  result->redim(qyqy, nen[i++], 1, 2);
  result->redim(qzqz, nsd == 3 ? nen[i++] : 0, 1, 3);

  result->redim(pp, nen[i], 1, -4);
  if (mat && mat->getPermeability())
  {
    result->redimOffDiag(qxqy, 0);
    result->redimOffDiag(qyqx, 0);
  }
  result->redimOffDiag(qxp, -1);
  result->redimOffDiag(qyp, -1);
  if (nsd == 3 && mat && mat->getPermeability())
  {
    result->redimOffDiag(qxqz, 0);
    result->redimOffDiag(qyqz, 0);
    result->redimOffDiag(qzqx, 0);
    result->redimOffDiag(qzqy, 0);
  }
  if (nsd == 3)
    result->redimOffDiag(qzp, -1);

  result->finalize();
  result->withLHS = true;
  result->rhsOnly = neumann || !calcMats || m_mode >= SIM::RHS_ONLY;
  return result;
}


bool CompatibleDarcy::initElement (const IntVec& MNPC,
                                   const UintVec& elem_sizes,
                                   const UintVec& basis_sizes,
                                   LocalIntegral& elmInt)
{
  return this->getElementSol(MNPC, elem_sizes, basis_sizes,
                             elmInt.vec, this->getNoSolutions());
}


bool CompatibleDarcy::getElementSol (const IntVec& MNPC,
                                     const UintVec& elem_sizes,
                                     const UintVec& basis_sizes,
                                     Vectors& eV, size_t nLevels) const
{
  int ierr = 0;
  if (primsol.empty() || primsol.front().empty())
    nLevels = 0;

  eV.resize(npv*nLevels);

  size_t j = 0;
  for (size_t level = 0; level < nLevels; ++level)
  {
    size_t ofs = 0;
    IntVec::const_iterator first = MNPC.begin();
    for (size_t i = 0; i < nsd; ++i)
    {
      IntVec::const_iterator last = first + elem_sizes[i];
      ierr += utl::gather(IntVec(first,last), 0,1, primsol[level],
                          eV[j++], ofs,ofs);
      first = last;
      ofs += basis_sizes[i];
    }

    IntVec MNPC2(first,first+elem_sizes[nsd]);
    ierr += utl::gather(MNPC2, 0,1, primsol[level],
                        eV[j++], ofs,ofs);
  }

  if (ierr > 0)
    std::cerr <<" *** CompatibleDarcy::getElementSol: Detected "
              << ierr <<" node numbers out of range."<< std::endl;

  return ierr == 0;
}


bool CompatibleDarcy::initElementBou (const IntVec& MNPC,
                                      const UintVec& elem_sizes,
                                      const UintVec& basis_sizes,
                                      LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool CompatibleDarcy::evalIntMx (LocalIntegral& elmInt,
                                 const MxFiniteElement& fe,
                                 const Vec3& X) const
{
  using IdxVec = std::array<int,3>;
  using IdxMat = std::array<IdxVec,3>;

  static constexpr IdxVec fluxDiagIdx{ qxqx, qyqy, qzqz };

  static constexpr IdxMat fluxMatIdx{{
      { qxqx, qxqy, qxqz },
      { qyqx, qyqy, qyqz },
      { qzqx, qzqy, qzqz }
    }};

  using CompatibleOps = CompatibleOperators::Weak;
  using ScalarOps     = EqualOrderOperators::Weak;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Matrix Kinv;
  const double mu = mat->getViscosity();
  const Vec3 K = mat->getPermeability(X);
  const bool hasKmat = mat->getPermeability(&Kinv);

  if (!elMat.A.empty() && calcMats)
  {
    if (hasKmat)
    {
      // The flux mass term is (mu K^-1 q, v), coupling the flux components
      if (Kinv.inverse() <= 0.0) return false;
      CompatibleOps::MassCoeff(elMat.A, Kinv, fe, fluxMatIdx, mu);
    }
    else for (size_t i = 0; i < nsd; ++i)
      ScalarOps::Mass(elMat.A[fluxDiagIdx[i]], fe, mu/K[i], i);

    CompatibleOps::Gradient(elMat.A, fe, { qxp, qyp, qzp });
  }

  Vec3 body(gravity);
  if (bodyforce)
    body += (*bodyforce)(X);
  if (!body.isZero())
  {
    body *= this->getDensity(fe);
    CompatibleOps::Source(elMat.b, fe, body, { Fqx, Fqy, Fqz });
  }

  if (source)
    ScalarOps::Source(elMat.b[Fp], fe, (*source)(X), 1, scalarBasis);

  return true;
}


bool CompatibleDarcy::evalBouMx (LocalIntegral& elmInt,
                                 const MxFiniteElement& fe,
                                 const Vec3& X, const Vec3& normal) const
{
  Vec3 fluxVec;
  if (flux)
    fluxVec = (*flux)(X)*normal;
  else if (vflux)
    fluxVec = (*vflux)(X);
  else if (tflux)
    fluxVec = (*tflux)(X,normal);
  else
  {
    std::cerr <<" *** CompatibleDarcy::evalBou: No fluxes."<< std::endl;
    return false;
  }

  CompatibleOperators::Weak::Source(static_cast<ElmMats&>(elmInt).b,
                                    fe, fluxVec, { Fqx, Fqy, Fqz });
  return true;
}


bool CompatibleDarcy::evalSol (Vector& s, const MxFiniteElement& fe,
                               const Vec3& X, const IntVec& MNPC,
                               const UintVec& elem_sizes,
                               const UintVec& basis_sizes) const
{
  s.clear();
  if (fluxOnly)
  {
    Vectors eV;
    if (!this->getElementSol(MNPC,elem_sizes,basis_sizes,eV))
      return false;

    s.resize(nsd);
    for (size_t i = 0; i < nsd; ++i)
      s[i] = eV[i].dot(fe.basis(i+1));
    return true;
  }
  else
    s.reserve(2+nsd);

  if (source)
    s.push_back((*source)(X));
  if (mat->isPorosityFunc())
    s.push_back(mat->getPorosity(X));

  if (mat->isPermeabilityFunc())
  {
    const Vec3 K = mat->getPermeability(X);
    s.push_back(K.ptr(), K.ptr()+nsd);
  }

  return true;
}


size_t CompatibleDarcy::getNoFields (int fld) const
{
  if (fld < 2)
    return npv;
  else if (fluxOnly)
    return nsd;

  if (size_t nfl = this->Darcy::getNoFields(fld); nfl > nsd)
    return nfl-nsd;

  return 0;
}


std::string CompatibleDarcy::getField1Name (size_t i, const char* pfx) const
{
  if (i == 11)
    return nsd == 2 ? "q_x&&q_y&&p" : "q_x&&q_y&&q_z&&p";
  if (i >= npv)
    return "";

  static const char* names[4] = {"q_x", "q_y", "q_z", "p"};
  const size_t index = i < nsd ? i : 3;
  return pfx ? pfx + std::string(" ") + names[index] : names[index];
}


std::string CompatibleDarcy::getField2Name (size_t i, const char* prefix) const
{
  return this->Darcy::getField2Name(fluxOnly ? i : nsd+i, prefix);
}


NormBase* CompatibleDarcy::getNormIntegrand (AnaSol* asol) const
{
  return new CompatibleDarcyNorm(*const_cast<CompatibleDarcy*>(this),asol);
}


CompatibleDarcyNorm::CompatibleDarcyNorm (CompatibleDarcy& p, AnaSol* a)
  : NormBase(p), flux(a ? a->getVectorSol() : nullptr),
    pressure(a ? a->getScalarSol() : nullptr)
{
  // Currently no secondary solution to recover for the mixed integrand,
  // and therefore no error estimates
}


bool CompatibleDarcyNorm::evalIntMx (LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const Vec3& X) const
{
  constexpr int H1_Ph   = DarcyNorm::H1_Ph;
  constexpr int H1_P    = DarcyNorm::H1_P;
  constexpr int H1_E_Ph = DarcyNorm::H1_E_Ph;
  constexpr int EXTENER = DarcyNorm::EXT_ENERGY;

  Darcy& problem = static_cast<Darcy&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const size_t nsd = problem.getNoSpaceDim();

  // Evaluate the inverse constitutive matrix at this point
  Matrix Kinv;
  if (!problem.getInvPermeability(X,Kinv))
    return false;

  // Evaluate the finite element Darcy velocity field and the pressure field
  Vector qh(nsd);
  for (size_t i = 0; i < nsd; ++i)
    qh[i] = pnorm.vec[i].dot(fe.basis(i+1));
  const double ph = pnorm.vec[nsd].dot(fe.basis(nsd+1));

  // Integrate the energy norm a(q^h,q^h)
  pnorm[H1_Qh] += qh.dot(Kinv*qh)*fe.detJxW;
  // Integrate the pressure L2-norm (p^h,p^h)
  pnorm[H1_Ph] += ph*ph*fe.detJxW;

  // Integrate the external energy (h,p^h)
  if (problem.extEner)
    pnorm[EXTENER] += problem.getPotential(X)*ph*fe.detJxW;

  if (flux)
  {
    // Evaluate the analytical Darcy velocity
    Vector q((*flux)(X).ptr(),nsd);
    // Integrate the energy norm a(q,q)
    pnorm[H1_Q] += q.dot(Kinv*q)*fe.detJxW;
    // Integrate the error in energy norm a(q-q^h,q-q^h)
    Vector error = q - qh;
    pnorm[H1_E_Qh] += error.dot(Kinv*error)*fe.detJxW;
  }

  if (pressure)
  {
    // Evaluate the analytical pressure
    const double p = (*pressure)(X);
    // Integrate the L2-norm (p,p)
    pnorm[H1_P] += p*p*fe.detJxW;
    // Integrate the pressure error in L2-norm (p-p^h,p-p^h)
    pnorm[H1_E_Ph] += (p-ph)*(p-ph)*fe.detJxW;
  }

  return true;
}


bool CompatibleDarcyNorm::evalBouMx (LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const Vec3& X, const Vec3& normal) const
{
  const Darcy& problem = static_cast<const Darcy&>(myProblem);
  if (!problem.extEner) return true;

  const size_t nsd = problem.getNoSpaceDim();
  const double flux = problem.getFlux(X,normal);
  const double ph = elmInt.vec[nsd].dot(fe.basis(nsd+1));

  // Integrate the external energy (h,p^h)
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  pnorm[DarcyNorm::EXT_ENERGY] += flux*ph*fe.detJxW;
  return true;
}


size_t CompatibleDarcyNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return group == 1 ? NUM_NORM : 0;
}


std::string CompatibleDarcyNorm::getName (size_t, size_t j,
                                          const char* prefix) const
{
  static const char* n[NUM_NORM] = {
    "(p^h,p^h)^0.5",
    "(h,p^h)^0.5",
    "(c^h,c^h)^0.5",
    "(p,p)^0.5",
    "(e,e)^0.5, e=p-p^h",
    "(c,c)^0.5",
    "(e,e)^0.5, e=c-c^h",
    "(e,e)^0.5, e=(p,c)-(p,c)^h",
    "a(q^h,q^h)^0.5",
    "a(q,q)^0.5",
    "a(e,e)^0.5, e=q-q^h"
  };

  if (!prefix)
    return n[j-1];

  return prefix + std::string(" ") + n[j-1];
}


bool CompatibleDarcyNorm::hasElementContributions (size_t i, size_t j) const
{
  if (i != 1)
    return false;

  if (j == 4 || j == 5)
    return pressure != nullptr;
  else if (j >= 2 && j <= 8)
    return false;
  else if (j == 10 || j == 11)
    return flux != nullptr;

  return true;
}
