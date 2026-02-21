// $Id$
//==============================================================================
//!
//! \file DarcyTransportCorr.C
//!
//! \date Jan 14 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy transport correction problem.
//!
//==============================================================================

#include "DarcyTransportCorr.h"

#include "AnaSol.h"
#include "BlockElmMats.h"
#include "ElmNorm.h"
#include "EqualOrderOperators.h"
#include "ExprFunctions.h"
#include "FiniteElement.h"
#include "Functions.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"

#include <array>
#include <cmath>
#include <cstring>


DarcyTransportCorr::DarcyTransportCorr (unsigned short int n, unsigned char n2)
  : IntegrandBase(n)
{
  npv = n;
  nf2 = n2;

  qq = 1;
  ql = n2 > 1 ? 4 : (n2 > 0 ? 3 : 0);
  qm = n2 > 1 ? 5 : 0;
  nM = n2 > 1 ? 10 : (n2 > 0 ? 5 : 3);

  Fq = 1;
  Fl = n2 > 0 ? 2 : 0;
  Fm = n2 > 1 ? 3 : 0;
  nV = n2 > 0 ? 2+n2 : 3;
}


DarcyTransportCorr::~DarcyTransportCorr() = default;


bool DarcyTransportCorr::parse (const tinyxml2::XMLElement* elem)
{
  const char* input = nullptr;
  if ((input = utl::getValue(elem,"observed_concentration")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tObserved concentration function: "<< input << std::endl;
    if (type == "expression")
      observed_C.reset(utl::parseExprRealFunc(input,true));
    else
      observed_C.reset(utl::parseRealFunc(input));
  }
  else if ((input = utl::getValue(elem,"input_source")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tInput source: "<< input << std::endl;
    if (type == "expression")
      input_source = std::make_unique<EvalFunction>(input);
    else
      input_source.reset(utl::parseRealFunc(input));
  }
  else if ((input = utl::getValue(elem,"input_velocity")))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tInput velocity function: "<< input << std::endl;
    if (type == "expression")
      input_q.reset(utl::parseExprVecFunc(input,true));
    else
      input_q = std::make_unique<VecFuncExpr>(input);
  }
  else if (!strcasecmp(elem->Value(),"penalty"))
  {
    utl::getAttribute(elem, "alpha", alpha);
    utl::getAttribute(elem, "beta", beta);
    utl::getAttribute(elem, "eps", eps);
  }
  else
    return false;

  return true;
}


LocalIntegral* DarcyTransportCorr::getLocalIntegral (size_t nen,
                                                     size_t, bool) const
{
  ElmMats* result = new ElmMats();

  result->resize(nM, nV);
  result->redim(nen*nsd);

  return result;
}


LocalIntegral*
DarcyTransportCorr::getLocalIntegral (const std::vector<size_t>& nen,
                                      size_t, bool) const
{
  BlockElmMats* result = new BlockElmMats(1+nf2,2);

  result->resize(nM, nV);
  result->redim(1, nen[0], nsd);
  for (size_t i = 0; i < nf2; i++)
    result->redim(2+i, nen[1], 1, -2);
  if (ql > 0)
    result->redimOffDiag(ql, -1);
  if (qm > 0)
    result->redimOffDiag(qm, -1);
  result->finalize();

  return result;
}


bool DarcyTransportCorr::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double  C   = (*observed_C)(X);
  const Vec3   dC   = observed_C->gradient(X);
  const double dCdt = observed_C->timeDerivative(X);

  const Vec3   q = (*input_q)(X);
  const double f = (*input_source)(X);

  const double scale = eps > 0.0 ? 1.0 / (eps + q.length2()) : 1.0;

  EqualOrderOperators::Weak::Mass(elMat.A[0], fe, scale);
  EqualOrderOperators::Weak::Source(elMat.b[0], fe, q, scale);

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (unsigned short int k = 1; k <= nsd; ++k)
        for (unsigned short int l = 1; l <= nsd; ++l)
        {
          const double mass_term = fe.dNdX(i,k) * fe.dNdX(j,l);
          const double transport_term = ((dC(k)*fe.N(i) + C*fe.dNdX(i,k)) *
                                         (dC(l)*fe.N(j) + C*fe.dNdX(j,l)));
          elMat.A[1]((i-1)*nsd+k, (j-1)*nsd+l) += mass_term * fe.detJxW;
          elMat.A[2]((i-1)*nsd+k, (j-1)*nsd+l) += transport_term * fe.detJxW;
        }

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (unsigned short int d = 1; d <= nsd; ++d)
    {
      const double transport_term = (f-dCdt) * (dC(d)*fe.N(i) + C*fe.dNdX(i,d));
      elMat.b[2]((i-1)*nsd + d) += transport_term * fe.detJxW;
    }

  return true;
}


bool DarcyTransportCorr::evalIntMx (LocalIntegral& elmInt,
                                    const MxFiniteElement& fe,
                                    const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double  C   = (*observed_C)(X);
  const Vec3   dC   = observed_C->gradient(X);
  const double dCdt = observed_C->timeDerivative(X);

  const Vec3   q = (*input_q)(X);
  const double R = (*input_source)(X);//this->residual(X);

  const double scale = 1.0;

  if (qq > 0)
    EqualOrderOperators::Weak::Mass(elMat.A[qq], fe, scale);

  if (ql > 0)
    for (size_t i = 1; i <= fe.basis(2).size(); ++i)
      for (size_t j = 1; j <= fe.basis(1).size(); ++j)
        for (unsigned short int d = 1; d <= nsd; ++d)
          elMat.A[ql]((j-1)*nsd+d,i) += fe.grad(1)(j,d) * fe.basis(2)(i) * fe.detJxW;

  if (qm > 0)
    for (size_t i = 1; i <= fe.basis(2).size(); ++i)
      for (size_t j = 1; j <= fe.basis(1).size(); ++j)
        for (unsigned short int d = 1; d <= nsd; ++d)
          elMat.A[qm]((j-1)*nsd+d,i) += (C*fe.grad(1)(j,d) + dC(d)*fe.basis(1)(j)) * fe.basis(2)(i) * fe.detJxW;

  if (Fq > 0)
    EqualOrderOperators::Weak::Source(elMat.b[Fq], fe, q);

  if (Fm > 0)
    EqualOrderOperators::Weak::Source(elMat.b[Fm], fe, R-dCdt, 1, 2);

  return true;
}


bool DarcyTransportCorr::finalizeElement (LocalIntegral& A)
{
  if (alpha > 0.0 && beta > 0.0 && ql == 0)
  {
    ElmMats& E = static_cast<ElmMats&>(A);

    E.A.front().add(E.A[1], alpha);
    E.A.front().add(E.A[2], beta);

    E.b.front().add(E.b[1], alpha);
    E.b.front().add(E.b[2], beta);
  }

  return true;
}


bool DarcyTransportCorr::evalSol2 (Vector& s, const Vectors& eV,
                                   const FiniteElement& fe, const Vec3& X) const
{
  s.resize(this->getNoFields(2));
  s[0] = (*observed_C)(X);
  s[1] = this->residual(X);
  s[2] = (*input_source)(X);
  s[3] = this->residual(eV.front(),fe.N,fe.dNdX,X);
  const Vec3 qh = (*input_q)(X);
  s[4] = (qh - this->evalSol(eV.front(),fe.N)).length();
  for (size_t i = 0; i < nsd; ++i)
    s[5+i] = qh[i];
  return true;
}


std::string DarcyTransportCorr::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    switch (nsd) {
    case 1: return "q_x";
    case 2: return "q_x&&q_y";
    case 3: return "q_x&&q_y&&q_z";
    }
  else if (i == 12)
    switch (nf2) {
    case 1: return "lambda";
    case 2: return "lambda&&mu";
    }

  std::string name;
  if (i < nsd)
    name = "q_" + std::string(1,'x'+i);
  else if (i-nsd == 0)
    name = "lambda";
  else if (i-nsd == 1)
    name = "mu";
  else
    name = "u" + std::to_string(i);

  return prefix ? prefix + std::string(" ") + name : name;
}


std::string DarcyTransportCorr::getField2Name (size_t i, const char* prefix) const
{
  const auto names3 = std::array {
    "observed_C",
    "residual",
    "source",
    "FE residual",
    "|q-q^|",
    "input q_x",
    "input q_y",
    "input q_z"
  };

  return prefix ? prefix + std::string(" ") + names3[i] : names3[i];
}


double DarcyTransportCorr::residual (const Vec3& X) const
{
  const double C = (*observed_C)(X);
  const Vec3 q = (*input_q)(X);
  const Vec3 grad_C = observed_C->gradient(X);
  const Tensor grad_q = input_q->gradient(X);
  const double f = (*input_source)(X);

  return f - (grad_q.trace()*C + q*grad_C);
}


double DarcyTransportCorr::residual (const Vector& eV,
                                     const Vector& N, const Matrix& dNdX,
                                     const Vec3& X) const
{
  const double  f   = (*input_source)(X);
  const double  C   = (*observed_C)(X);
  const Vec3   dCdX = observed_C->gradient(X);
  const Vec3    q   = this->evalSol(eV,N);
  const Tensor dqdX = this->evalGrd(eV,dNdX);

  return f - (dqdX.trace()*C + q*dCdX);
}


Vec3 DarcyTransportCorr::evalSol (const Vector& eV, const Vector& N) const
{
  Vec3 q;
  for (unsigned short int i = 0; i < nsd; i++)
    q[i] = eV.dot(N,i,nsd);

  return q;
}


Tensor DarcyTransportCorr::evalGrd (const Vector& eV, const Matrix& dNdX) const
{
  Tensor dqdX(nsd);
  for (unsigned short int j = 1; j <= nsd; j++)
    for (unsigned short int i = 1; i <= nsd; i++)
      dqdX(i,j) = eV.dot(dNdX.ptr(j-1),dNdX.rows(),i-1,nsd);

  return dqdX;
}


Vec3 DarcyTransportCorr::evalInput (const Vec3& X) const
{
  return (*input_q)(X);
}


double DarcyTransportCorr::divInput (const Vec3& X) const
{
  return input_q->gradient(X).trace();
}


double DarcyTransportCorr::evalTracer (const Vec3& X) const
{
  return (*observed_C)(X);
}


NormBase* DarcyTransportCorr::getNormIntegrand (AnaSol* asol) const
{
  return new DarcyTCNorm(*const_cast<DarcyTransportCorr*>(this));
}


bool DarcyTCNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                           const TimeDomain& time, const Vec3& X) const
{
  DarcyTransportCorr& problem = static_cast<DarcyTransportCorr&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  const double c_i = problem.evalTracer(X);
  const Vec3 q_inp = problem.evalInput(X);
  const Vec3 q_hat = problem.evalSol(pnorm.vec.front(),fe.N);
  const Tensor grd = problem.evalGrd(pnorm.vec.front(),fe.dNdX);
  const double div = grd.trace();
  const double res = problem.residual(pnorm.vec.front(),fe.N,fe.dNdX,X);

  pnorm[0] += q_hat.length2() * fe.detJxW;
  pnorm[1] += c_i*q_hat.length2() * fe.detJxW;
  pnorm[2] += (q_hat-q_inp).length2() * fe.detJxW;
  pnorm[3] += div*div * fe.detJxW;
  pnorm[4] += res*res * fe.detJxW;

  const double divQ = problem.divInput(X);
  pnorm[5] += q_inp.length2() * fe.detJxW;
  pnorm[6] += divQ*divQ * fe.detJxW;

  return true;
}


bool DarcyTCNorm::evalIntMx (LocalIntegral& elmInt, const MxFiniteElement& fe,
                             const TimeDomain& time, const Vec3& X) const
{
  return this->evalInt(elmInt,fe,time,X);
}


size_t DarcyTCNorm::getNoFields (int group) const
{
  return group < 1 ? 1 : (group > 1 ? 0 : 7);
}


std::string DarcyTCNorm::getName (size_t, size_t j, const char* prefix) const
{
  const auto names = std::array {
    "L2(q^_h)",
    "L2(c*q^_h)",
    "L2(q^_h-q_h)",
    "L2(div q^h)",
    "L2(residual)",
    "L2(q)",
    "L2(div q)"
  };

  return prefix ? prefix + std::string(" ") + names[j-1] : names[j-1];
}
