// $Id$
//==============================================================================
//!
//! \file MixedDarcy.C
//!
//! \date oct 20 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mixed Darcy flow problems.
//!
//==============================================================================

#include "MixedDarcy.h"

#include "AnaSol.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "BlockElmMats.h"
#include "ElmNorm.h"
#include "Fields.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "SIMbase.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <cmath>
#include <ext/alloc_traits.h>
#include <iostream>
#include <memory>
#include <vector>


MixedDarcy::MixedDarcy (unsigned short int n, int torder) :
  Darcy(n, torder)
{
  pp = 1;
  cc = 2;
  cp = 4;
  sourceC = nullptr;
  if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2 ||
      ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS2)
    cBasis = 2;
  else {
    npv = 2;
    cBasis = 1;
  }
}


LocalIntegral* MixedDarcy::getLocalIntegral (size_t nen, size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2, 1);

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(5, 3);
  result->redim(pp, nen, 1, 1);
  result->redim(cc, nen, 1, 1);
  result->redimOffDiag(cp, 0);
  result->redimNewtonMat();

  return result;
}


LocalIntegral* MixedDarcy::getLocalIntegral (const std::vector<size_t>& nen,
                                             size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2, 2);

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(5, 3);
  result->redim(pp, nen[0], 1, 1);
  result->redim(cc, nen[1], 1, 2);
  result->redimOffDiag(cp, 0);
  result->redimNewtonMat();

  return result;
}


bool MixedDarcy::initElement (const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes,
                              const std::vector<size_t>& basis_sizes,
                              LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(2*this->getNoSolutions());
  int ierr = 0;
  auto fstart = MNPC.begin() + elem_sizes[0];
  auto fend = fstart + elem_sizes[1];
  for (size_t k = 0; k < elmInt.vec.size() / 2; ++k) {
    ierr  += utl::gather(IntVec(MNPC.begin(), fstart), 1, primsol[k], elmInt.vec[k*2+0]);
    ierr +=  utl::gather(IntVec(fstart, fend), 0, 1, primsol[k],
                         elmInt.vec[k*2+1], basis_sizes[0], basis_sizes[0]);
  }

  if (ierr != 0)
    std::cerr << " *** MixedDarcy::initElement: Detected " << ierr/2
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool MixedDarcy::initElement (const std::vector<int>& MNPC,
                              LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(2*this->getNoSolutions());
  int ierr = 0;
  for (size_t k = 0; k < elmInt.vec.size() / 2; ++k) {
    Matrix tmp(2,MNPC.size());
    utl::gather(MNPC,npv,primsol[k],tmp);
    elmInt.vec[2*k] = tmp.getRow(1);
    elmInt.vec[2*k+1] = tmp.getRow(2);
  }

  if (ierr != 0)
    std::cerr << " *** MixedDarcy::initElement: Detected " << ierr
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool MixedDarcy::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                          const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!this->Darcy::evalInt(elmInt, fe, time, X))
    return false;

  if (!dispersivity || !porosity)
    return false;

  const double D = (*this->dispersivity)(X);
  const double phi = (*this->porosity)(X);

  WeakOps::Laplacian(elMat.A[cc], fe, D, false, cBasis);

  double cn = this->concentration(elmInt.vec, fe, 0);

  Vec3 perm = this->getPermeability(X);
  elMat.A[cp].multiply(fe.grad(cBasis),fe.grad(1),false,true,true,perm[0]*cn*fe.detJxW);

  if (sourceC)
    WeakOps::Source(elMat.b[2], fe, (*sourceC)(X), 1, cBasis);

  if (bdf.getActualOrder() > 0) {
    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = this->concentration(elmInt.vec, fe, t);
      c -= val * phi * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b[2], fe, c, 1, cBasis);
    WeakOps::Mass(elMat.A[cc], fe, phi*bdf[0] / time.dt, cBasis);
  }

  return true;
}


bool MixedDarcy::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                          const Vec3& X, const Vec3& normal) const
{
  if (!this->Darcy::evalBou(elmInt,fe,X,normal))
    return false;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!dispersivity)
    return false;

  const double D = (*this->dispersivity)(X);

  Vec3 perm = this->getPermeability(X);

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int k = 1; k <= nsd; ++k)
        elMat.A[cp](i,j) += (fe.N(j)*perm[k-1]*fe.dNdX(j,k)*normal[k-1] - D*fe.dNdX(j,k))*fe.N(i)*fe.detJxW;

  return true;
}


bool MixedDarcy::finalizeElement (LocalIntegral& A)
{
  if (m_mode == SIM::RHS_ONLY && !reacInt) {
    ElmMats& elMat = static_cast<ElmMats&>(A);
    elMat.A[pp].multiply(elMat.vec[0], elMat.b[1], false, -1);
    elMat.A[cc].multiply(elMat.vec[1], elMat.b[2], false, -1);
    elMat.A[cp].multiply(elMat.vec[0], elMat.b[2], false, -1);
  }

  return true;
}


bool MixedDarcy::evalSolInt (Vector& s, const Vectors& eV,
                             const FiniteElement& fe, const Vec3& X) const
{
  if (!this->Darcy::evalSol2(s,eV,fe,X))
    return false;

  // Evaluate the concentration gradient
  Vec3 dCh = this->concentrationGradient(eV,fe,0);

  for (int i = 0; i < nsd; ++i)
    s.push_back(dCh[i]);

  return true;
}


bool MixedDarcy::evalSol (Vector& s, const FiniteElement& fe,
                          const Vec3& X,
                          const std::vector<int>& MNPC) const
{
  ElmMats A;
  if (!const_cast<MixedDarcy*>(this)->initElement(MNPC,A))
    return false;

  return this->evalSolInt(s,A.vec,fe,X);
}

bool MixedDarcy::evalSol (Vector& s,
                          const MxFiniteElement& fe,
                          const Vec3& X,
                          const std::vector<int>& MNPC,
                          const std::vector<size_t>& elem_sizes,
                          const std::vector<size_t>& basis_sizes) const
{
  ElmMats A;
  if (!const_cast<MixedDarcy*>(this)->initElement(MNPC,elem_sizes,basis_sizes,A))
    return false;

  return this->evalSolInt(s,A.vec,fe,X);
}


std::string MixedDarcy::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return cBasis == 1 ? "p&&c" : "p";
  if (i == 12)
    return "c";

  if (i >= 2)
    return "";

  static const char* s[2] = {"p", "c"};
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string MixedDarcy::getField2Name (size_t i, const char* prefix) const
{
  if (nsd == 2 && i > 1)
    ++i;

  if (i >= 6) return "";

  static const char* s[6] = {"p,x","p,y","p,z",
                             "c,x","c,y","c,z"};

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


double MixedDarcy::concentration (const Vectors& vec,
                                  const FiniteElement& fe,
                                  size_t level) const
{
  return fe.basis(cBasis).dot(vec[level*2+1]);
}


Vec3 MixedDarcy::concentrationGradient (const Vectors& vec,
                                        const FiniteElement& fe,
                                        size_t level) const
{
  Vector dCh(nsd);
  fe.grad(cBasis).multiply(vec[2*level+1],dCh,true);

  return dCh;
}


double MixedDarcy::pressure (const Vectors& eV,
                             const FiniteElement& fe,
                             size_t level) const
{
  return fe.basis(1).dot(eV[2*level]);
}


Vec3 MixedDarcy::pressureGradient (const Vectors& eV,
                                   const FiniteElement& fe,
                                   size_t level) const
{
  Vector dPh(nsd);
  fe.grad(1).multiply(eV[2*level],dPh,true);

  return dPh;
}


void MixedDarcy::getSolutionNorms (const SIMbase& sim,
                                   const Vector& solution,
                                   double& dNorm,
                                   double* dMax, size_t* iMax) const
{
  if (ASMmxBase::Type == ASMmxBase::REDUCED_CONT_RAISE_BASIS2 ||
      ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS2) {
    double dNormp = sim.solutionNorms(solution, dMax, iMax, 1, 'D');
    double dNormc = sim.solutionNorms(solution, dMax+1, iMax+1, 1, 'P');

    dNorm = sqrt(pow(dNormp,2.0)+pow(dNormc,2.0));
  } else
    dNorm = sim.solutionNorms(solution,dMax,iMax,2);
}



NormBase* MixedDarcy::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new MixedDarcyNorm(*const_cast<MixedDarcy*>(this),
                              asol->getScalarSecSol(0),
                              asol->getScalarSecSol(1));

  else
    return new MixedDarcyNorm(*const_cast<MixedDarcy*>(this));
}


MixedDarcyNorm::MixedDarcyNorm (MixedDarcy& p, VecFunc* a, VecFunc* c)
  : DarcyNorm(p,a), anac(c)
{
}


bool MixedDarcyNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  if (!this->DarcyNorm::evalInt(elmInt,fe,X))
    return false;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const MixedDarcy& problem = static_cast<const MixedDarcy&>(myProblem);
  const double D = problem.getDispersivity(X);

  // Evaluate the concentration field gradient
  Vec3 dCh = problem.concentrationGradient(elmInt.vec, fe, 0);

  pnorm[H1_Ch] += D*dCh*dCh*fe.detJxW;

  Vec3 dC;
  if (anac) {
    dC = (*anac)(X);
    pnorm[H1_C] += D*dC*dC*fe.detJxW;
    Vec3 error =  dC-dCh;
    double E = D*error*error*fe.detJxW;
    pnorm[H1_E_Ch] += E;
    pnorm[TOTAL_NORM_E] += E;
  }

  size_t ip = this->getNoFields(1);
  size_t f = 2;
  for (const Vector& psol : pnorm.psol) {
    if (!projFields.empty() || !psol.empty())
    {
      Vec3 dCr;
      if (!projFields.empty() && projFields[f-2]) {
        Vector vals;
        projFields[f-2]->valueFE(fe, vals);
        for (size_t i = 0; i < fe.dNdX.cols(); ++i)
          dCr[i] = vals[i+fe.dNdX.cols()];
      } else {
        // Evaluate projected concentration field
        for (size_t j = 0; j < fe.dNdX.cols(); j++)
          dCr[j] = psol.dot(fe.N,j+fe.dNdX.cols(),nrcmp);
      }

      // Integrate the energy norm a(c^r,c^r)
      pnorm[ip+H1_Cr] += D*dCr*dCr*fe.detJxW;
      // Integrate the estimated error in energy norm a(c^r-c^h,c^r-c^h)
      Vec3 error = dCr - dCh;
      double E = D*error*error*fe.detJxW;
      pnorm[ip+H1_Cr_Ch] += E;
      pnorm[ip+TOTAL_NORM_REC] += E;

      if (anac)
      {
        // Integrate the error in the projected solution a(c-c^r,c-c^r)
        error = dC - dCr;
        E = D*error*error*fe.detJxW;
        pnorm[ip+H1_E_Cr] += E;
        pnorm[ip+TOTAL_E_REC] += E;
      }
    }
    ip += this->getNoFields(f++);
  }

  return true;
}
