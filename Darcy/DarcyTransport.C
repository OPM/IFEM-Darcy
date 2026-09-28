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
#include "DarcyMaterial.h"

#include "AnaSol.h"
#include "BlockElmMats.h"
#include "ElmNorm.h"
#include "EqualOrderOperators.h"
#include "Fields.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"


DarcyTransport::DarcyTransport (unsigned short int n, int torder) :
  Darcy(n, torder)
{
  pp = 1;
  cc = 2;
  cp = 4;
  npv = 2;
}


DarcyTransport::~DarcyTransport() = default;


LocalIntegral* DarcyTransport::getLocalIntegral (size_t nen, size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2, 1);

  result->rhsOnly = neumann || !calcMats;
  result->withLHS = !neumann;
  result->resize(5, 3);
  result->redim(pp, nen, 1, 1);
  result->redim(cc, nen, 1, 1);
  result->redimOffDiag(cp, 0);
  result->finalize();

  return result;
}


void DarcyTransport::setCSource (RealFunc* s)
{
  sourceC.reset(s);
}


bool DarcyTransport::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  return this->getElementSol(MNPC,elmInt.vec,this->getNoSolutions());
}


bool DarcyTransport::getElementSol (const std::vector<int>& MNPC,
                                    Vectors& eV, size_t nSol) const
{
  int ierr = 0;
  if (primsol.empty() || primsol.front().empty())
    nSol = 0;

  // Extract the element-level solution vectors.
  // Here we put each solution component as separate element-level vectors.
  eV.clear();
  eV.reserve(npv*nSol);
  for (size_t k = 0; k < nSol && !primsol[k].empty(); ++k)
  {
    Matrix tmp;
    ierr += utl::gather(MNPC,npv,primsol[k],tmp);
    for (unsigned short int i = 1; i <= npv; i++)
      eV.push_back(tmp.getRow(i));
  }

  if (ierr > 1)
    std::cerr <<" *** DarcyTransport::getElementSol: Detected "<< ierr/2
              <<" node numbers out of range."<< std::endl;

  return ierr == 0;
}


bool DarcyTransport::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  if (!this->Darcy::evalInt(elmInt, fe, time, X))
    return false;

  using WeakOps = EqualOrderOperators::Weak;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double D = mat->getDispersivity(X);
  const double phiDt = mat->getPorosity(X) / time.dt;

  if (!elMat.A.empty() && calcMats) {
    WeakOps::Laplacian(elMat.A[cc], fe, D, false);

    Matrix Kmat;
    double cn = this->concentration(elmInt.vec, fe, 0);
    if (mat->getPermeability(&Kmat))
      cn *= Kmat(1,1);
    else
      cn *= mat->getPermeability(X).x;
    elMat.A[cp].multiply(fe.dNdX,fe.dNdX,false,true,true,cn*fe.detJxW);
  }

  if (sourceC)
    WeakOps::Source(elMat.b[2], fe, (*sourceC)(X), 1);

  if (bdf.getActualOrder() > 0) {
    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++)
      c -= this->concentration(elmInt.vec,fe,t) * phiDt*bdf[t];
    WeakOps::Source(elMat.b[2], fe, c, 1);
    if (!elMat.A.empty() && calcMats)
      WeakOps::Mass(elMat.A[cc], fe, phiDt*bdf[0]);
  }

  return true;
}


bool DarcyTransport::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3& X, const Vec3& normal) const
{
  if (!this->Darcy::evalBou(elmInt,fe,X,normal) || !mat)
    return false;
  else if (mat->getPermeability())
  {
    std::cerr <<" *** DarcyTransport::evalBou():"
              <<" Implemented for diagonal permeability only."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double D = mat->getDispersivity(X);
  const Vec3 K = mat->getPermeability(X);

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int k = 1; k <= nsd; ++k)
        elMat.A[cp](i,j) += (fe.N(j)*K(k)*fe.dNdX(j,k)*normal(k) - D*fe.dNdX(j,k))*fe.N(i)*fe.detJxW;

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
  s.clear();
  s.reserve(3+3*nsd);

  Vectors eV;
  if (!this->getElementSol(MNPC,eV))
    return false;

  if (!this->evalDarcyVel(s,eV,fe,X))
    return false;
  else if (fluxOnly)
    return true;

  Vec3 dC = this->concentrationGradient(eV,fe,0);
  s.push_back(dC.ptr(),dC.ptr()+nsd);

  if (source)
    s.push_back((*source)(X));
  if (sourceC)
    s.push_back((*sourceC)(X));
  if (mat->isPorosityFunc())
    s.push_back(mat->getPorosity(X));

  if (mat->isPermeabilityFunc())
  {
    Vec3 K = mat->getPermeability(X);
    s.push_back(K.ptr(),K.ptr()+nsd);
  }

  return true;
}


size_t DarcyTransport::getNoFields (int fld) const
{
  if (fld < 2)
    return 2;
  else if (!mat)
    return 0;
  else if (fluxOnly)
    return nsd;

  size_t nField = nsd*(mat->isPermeabilityFunc() ? 3 : 2);
  if (source)
    ++nField;
  if (sourceC)
    ++nField;
  if (mat->isPorosityFunc())
    ++nField;

  return nField;
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
  if (i >= this->DarcyTransport::getNoFields(2)) return "";

  if (nsd == 2 && i > 1)
    ++i;
  if (nsd == 2 && i > 4)
    ++i;
  if (!source && i > 5)
    ++i;
  if (!sourceC && i > 6)
    ++i;
  if (!mat->isPorosityFunc() && i > 7)
    ++i;

  static const char* s[12] = { "v_x", "v_y", "v_z",
                               "c,x", "c,y", "c,z",
                               "source", "source_c",
                               "porosity",
                               "perm_x", "perm_y", "perm_z" };

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


double DarcyTransport::concentration (const Vectors& vec,
                                      const FiniteElement& fe,
                                      size_t level) const
{
  return fe.N.dot(vec[level*2+1]);
}


Vec3 DarcyTransport::concentrationGradient (const Vectors& vec,
                                            const FiniteElement& fe,
                                            size_t level) const
{
  RealArray dCh(nsd);
  fe.dNdX.multiply(vec[2*level+1],dCh,true);

  return dCh;
}


double DarcyTransport::pressure (const Vectors& eV,
                                 const FiniteElement& fe,
                                 size_t level) const
{
  return fe.N.dot(eV[2*level]);
}


NormBase* DarcyTransport::getNormIntegrand (AnaSol* asol) const
{
  return new DarcyTransportNorm(*const_cast<DarcyTransport*>(this),
                                asol ? asol->getScalarSecSol(0) : nullptr,
                                asol ? asol->getScalarSecSol(1) : nullptr);
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
  const double D = problem.getDispersivity(X);

  // Evaluate the concentration field gradient
  Vec3 dCh = problem.concentrationGradient(elmInt.vec, fe, 0);

  pnorm[H1_Ch] += D*dCh*dCh*fe.detJxW;

  double E = 0.0;
  Vec3 dC, dCr, error;
  if (anac)
  {
    dC = (*anac)(X);
    pnorm[H1_C] += D*dC*dC*fe.detJxW;
    error = dC - dCh;
    E = D*error*error*fe.detJxW;
    pnorm[H1_E_Ch] += E;
    pnorm[TOTAL_NORM_E] += E;
  }

  const size_t nsd = fe.dNdX.cols();
  size_t ip = this->getNoFields(1);
  for (size_t k = 0; k < pnorm.psol.size(); ip += this->getNoFields(2+k++))
    if (!prjFld.empty() || !pnorm.psol[k].empty())
    {
      // Evaluate projected concentration field
      if (prjFld.size() > k && prjFld[k])
      {
        Vector vals;
        prjFld[k]->valueFE(fe,vals);
        dCr = Vec3(vals.ptr()+nsd,nsd);
      }
      else
        for (size_t i = 0; i < nsd; i++)
          dCr[i] = pnorm.psol[k].dot(fe.N,nsd+i,nrcmp);

      // Integrate the energy norm a(c^r,c^r)
      pnorm[ip+H1_Cr] += D*dCr*dCr*fe.detJxW;
      // Integrate the estimated error in energy norm a(c^r-c^h,c^r-c^h)
      error = dCr - dCh;
      E = D*error*error*fe.detJxW;
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

  return true;
}
