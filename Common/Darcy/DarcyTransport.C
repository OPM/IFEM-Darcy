// $Id$
//==============================================================================
//!
//! \file DarcyTransport.C
//!
//! \date oct 20 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mixed Darcy flow problems.
//!
//==============================================================================

#include "DarcyTransport.h"

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


DarcyTransport::DarcyTransport (unsigned short int n, int torder) :
  Darcy(n, torder)
{
  pp = 1;
  cc = 2;
  cp = 4;
  sourceC = nullptr;
  npv = 2;
}


LocalIntegral* DarcyTransport::getLocalIntegral (size_t nen, size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2, 1);

  result->rhsOnly = neumann || this->reuseMats;
  result->withLHS = !neumann;
  result->resize(5, 3);
  result->redim(pp, nen, 1, 1);
  result->redim(cc, nen, 1, 1);
  result->redimOffDiag(cp, 0);
  result->finalize();

  return result;
}


bool DarcyTransport::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(2*this->getNoSolutions());
  int ierr = 0;
  for (size_t k = 0; k < elmInt.vec.size() / 2 && !primsol[k].empty(); ++k) {
    Matrix tmp(2,MNPC.size());
    utl::gather(MNPC,npv,primsol[k],tmp);
    elmInt.vec[2*k] = tmp.getRow(1);
    elmInt.vec[2*k+1] = tmp.getRow(2);
  }

  if (ierr != 0)
    std::cerr << " *** DarcyTransport::initElement: Detected " << ierr
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool DarcyTransport::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!this->Darcy::evalInt(elmInt, fe, time, X))
    return false;

  if (!this->mat)
    return false;

  const double D = this->getMaterial().getDispersivity(X);
  const double phi = this->getMaterial().getPorosity(X);

  if (!elMat.A.empty() && !reuseMats) {
    WeakOps::Laplacian(elMat.A[cc], fe, D, false);

    double cn = this->concentration(elmInt.vec, fe, 0);
    Vec3 perm = this->getMaterial().getPermeability(X);
    elMat.A[cp].multiply(fe.grad(1),fe.grad(1),false,true,true,perm[0]*cn*fe.detJxW);
  }

  if (sourceC)
    WeakOps::Source(elMat.b[2], fe, (*sourceC)(X), 1);

  if (bdf.getActualOrder() > 0) {
    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = this->concentration(elmInt.vec, fe, t);
      c -= val * phi * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b[2], fe, c, 1);
    if (!elMat.A.empty() && !reuseMats)
      WeakOps::Mass(elMat.A[cc], fe, phi*bdf[0] / time.dt);
  }

  return true;
}


bool DarcyTransport::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3& X, const Vec3& normal) const
{
  if (!this->Darcy::evalBou(elmInt,fe,X,normal))
    return false;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!this->mat)
    return false;

  const double D = this->getMaterial().getDispersivity(X);

  Vec3 perm = this->getMaterial().getPermeability(X);

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int k = 1; k <= nsd; ++k)
        elMat.A[cp](i,j) += (fe.N(j)*perm[k-1]*fe.dNdX(j,k)*normal[k-1] - D*fe.dNdX(j,k))*fe.N(i)*fe.detJxW;

  return true;
}


bool DarcyTransport::finalizeElement (LocalIntegral& A)
{
  if (m_mode == SIM::RHS_ONLY && !reacInt) {
    ElmMats& elMat = static_cast<ElmMats&>(A);
    elMat.A[pp].multiply(elMat.vec[0], elMat.b[1], false, -1);
    elMat.A[cc].multiply(elMat.vec[1], elMat.b[2], false, -1);
    elMat.A[cp].multiply(elMat.vec[0], elMat.b[2], false, -1);
  }

  return true;
}


bool DarcyTransport::evalSol (Vector& s, const FiniteElement& fe,
                              const Vec3& X,
                              const std::vector<int>& MNPC) const
{
  ElmMats A;
  if (!const_cast<DarcyTransport*>(this)->initElement(MNPC,A))
    return false;

  if (!this->Darcy::evalDarcyVel(s,A.vec,fe,X))
    return false;

  s.push_back(source ? (*source)(X) : 0.0);
  s.push_back(sourceC ? (*sourceC)(X) : 0.0);
  s.push_back(this->getMaterial().getPorosity(X));

  Vec3 dCh = this->concentrationGradient(A.vec,fe,0);
  for (int i = 0; i < nsd; ++i)
    s.push_back(dCh[i]);

  Vec3 perm = this->getMaterial().getPermeability(X);
  for (int i = 0; i < nsd; ++i)
    s.push_back(perm[i]);

  return true;
}

std::string DarcyTransport::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return "p&&c";

  if (i >= 2)
    return "";

  static const char* s[2] = {"p", "c"};
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string DarcyTransport::getField2Name (size_t i, const char* prefix) const
{
  if (i >= (nsd == 2 ? 9 : 12)) return "";

  static const char* s2[9] = {"v_x", "v_y",
                              "source", "source_c",
                              "porosity",
                              "c,x", "c,y",
                              "perm_x", "perm_y"};

  static const char* s3[12] = {"v_x","v_y","v_z",
                               "source","source_c",
                               "porosity",
                               "c,x","c,y","c,z",
                               "perm_x", "perm_y", "perm_z"};

  const char** s = (nsd == 2 ? s2 : s3);

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


double DarcyTransport::concentration (const Vectors& vec,
                                      const FiniteElement& fe,
                                      size_t level) const
{
  return fe.basis(1).dot(vec[level*2+1]);
}


Vec3 DarcyTransport::concentrationGradient (const Vectors& vec,
                                            const FiniteElement& fe,
                                            size_t level) const
{
  Vector dCh(nsd);
  fe.grad(1).multiply(vec[2*level+1],dCh,true);

  return dCh;
}


double DarcyTransport::pressure (const Vectors& eV,
                                 const FiniteElement& fe,
                                 size_t level) const
{
  return fe.basis(1).dot(eV[2*level]);
}


Vec3 DarcyTransport::pressureGradient (const Vectors& eV,
                                       const FiniteElement& fe,
                                       size_t level) const
{
  Vector dPh(nsd);
  fe.grad(1).multiply(eV[2*level],dPh,true);

  return dPh;
}



NormBase* DarcyTransport::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new DarcyTransportNorm(*const_cast<DarcyTransport*>(this),
                                  asol->getScalarSecSol(0),
                                  asol->getScalarSecSol(1));

  else
    return new DarcyTransportNorm(*const_cast<DarcyTransport*>(this));
}


DarcyTransportNorm::DarcyTransportNorm (DarcyTransport& p, VecFunc* a, VecFunc* c)
  : DarcyNorm(p,a), anac(c)
{
}


bool DarcyTransportNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                  const TimeDomain& time, const Vec3& X) const
{
  if (!this->DarcyNorm::evalInt(elmInt,fe,X))
    return false;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const DarcyTransport& problem = static_cast<const DarcyTransport&>(myProblem);
  const double D = problem.getMaterial().getDispersivity(X);

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
          dCr[i] = vals[i+3+fe.dNdX.cols()];
      } else {
        // Evaluate projected concentration field
        for (size_t i = 0; i < fe.dNdX.cols(); ++i)
          dCr[i] = psol.dot(fe.N,i+3+fe.dNdX.cols(),nrcmp);
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
