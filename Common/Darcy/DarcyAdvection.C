// $Id$
//==============================================================================
//!
//! \file DarcyAdvection.C
//!
//! \date Aug 26 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Darcy advection problems.
//!
//==============================================================================

#include "DarcyAdvection.h"
#include "DarcySolutions.h"
#include "Darcy.h"

#include "ElmMats.h"
#include "ElmNorm.h"
#include "ExprFunctions.h"
#include "Field.h"
#include "FiniteElement.h"
#include "IFEM.h"
#include "SIMbase.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"

#include <ext/alloc_traits.h>
#include "tinyxml2.h"


DarcyAdvection::DarcyAdvection (unsigned short int n,
                                const Darcy& drc1, int torder) :
  IntegrandBase(n), bdf(torder), drc(drc1)
{
  npv = 1;
  this->registerVector("pressure",&pVec);
  primsol.resize(1+torder);
}


bool DarcyAdvection::parse (const tinyxml2::XMLElement* elem)
{
  if (const char* input = utl::getValue(elem,"source"); input)
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tSource function:";

    RealFunc* src = nullptr;
    if (type == "expression")
    {
      IFEM::cout <<" "<< input << std::endl;
      src = new EvalFunction(input);
    }
    else if (type == "diracsum")
    {
      double tol = 1e-2;
      utl::getAttribute(elem,"pointTol",tol);
      if (input)
      {
        IFEM::cout << " DiracSum";
        DiracSum* f = new DiracSum(tol,nsd);
        if (f->parse(input))
          src = f;
        else
          delete f;
      }
    }
    else if (type == "elementsum" && ownerSim->createFEMmodel('y'))
    {
      IFEM::cout << " ElementSum";
      ElementSum* f = new ElementSum(nsd);
      if (f->parse(input,*ownerSim))
        src = f;
      else
        delete f;
    }

    if (src)
      source.reset(src);
  }
  else
    return false;

  return true;
}


LocalIntegral* DarcyAdvection::getLocalIntegral (size_t nen,
                                                 size_t iEl, bool neumann) const
{
  LocalIntegral* res = this->IntegrandBase::getLocalIntegral(nen,iEl,neumann);
  if (useLCache)
    static_cast<ElmMats*>(res)->rhsOnly = reuseMats;

  return res;
}


bool DarcyAdvection::initElement (const std::vector<int>& MNPC,
                                  const FiniteElement& fe,
                                  const Vec3& XC,
                                  size_t nPt, LocalIntegral& elmInt)
{
  if (fe.iel > 0 && reuseMats)
  {
    size_t iel = fe.iel - 1;
    ElmMats* A = dynamic_cast<ElmMats*>(&elmInt);
    if (A && iel < this->myKmats.size() && !A->A.empty())
      A->A[0] = this->myKmats[iel];
  }

  return this->IntegrandBase::initElement(MNPC,fe,XC,nPt,elmInt);
}


bool DarcyAdvection::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const double D = mat->getDispersivity(X);
  const double phi = mat->getPorosity(X);

  if (!elMat.A.empty() && !reuseMats) {
    WeakOps::Laplacian(elMat.A[0], fe, D, false);

    RealArray q;
    if (!this->evalDarcyVel(q,fe,X))
      return false;

    WeakOps::Advection(elMat.A[0], fe, q, 1.0, WeakOperators::CONSERVATIVE);
  }

  if (source)
    WeakOps::Source(elMat.b[0], fe, (*source)(X), 1);

  if (bdf.getActualOrder() > 0) {
    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = elmInt.vec[t].dot(fe.N);
      c -= val * phi * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b[0], fe, c, 1);
    if (!elMat.A.empty() && !reuseMats)
      WeakOps::Mass(elMat.A[0], fe, phi*bdf[0] / time.dt);
  }

  return true;
}


bool DarcyAdvection::finalizeElement (LocalIntegral& elmInt,
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


bool DarcyAdvection::evalDarcyVel (RealArray& q,
                                   const FiniteElement& fe, const Vec3& X) const
{
  Matrix K;
  drc.formKmatrix(K, X);

  Vector dP;
  this->pField->gradFE(fe, dP);

  return K.multiply(dP, q, - 1.0 / mat->getViscosity());
}


std::string DarcyAdvection::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return "c";

  if (!prefix) return "c";

  return prefix + std::string(" c");
}


std::string DarcyAdvection::getField2Name (size_t i, const char* prefix) const
{
  if (i >= this->getNoFields(2)) return "";

  static const char* s[4] = {"source_c", "c,x", "c,y", "c,z"};

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


bool DarcyAdvection::evalSol2 (Vector& s, const Vectors& eV,
                               const FiniteElement& fe, const Vec3& X) const
{
  s.clear();
  s.reserve(1+nsd);
  s.push_back(source ? (*source)(X) : 0.0);

  Vector dCh(nsd);
  fe.grad(1).multiply(eV.front(),dCh,true);
  for (int i = 0; i < nsd; ++i)
    s.push_back(dCh[i]);

  return true;
}


void DarcyAdvection::setNamedField (const std::string& name, Field* field)
{
  if (name == "pressure")
    pField.reset(field);
  else
    delete field;
}


void DarcyAdvection::initLHSbuffers (size_t nEl)
{
  if (!useLCache)
    return;

  if (nEl > 1) {
    this->myKmats.resize(nEl);
    this->reuseMats = false;
  } else if (nEl == 1)
    this->reuseMats = false;
  else if (!myKmats.empty())
    this->reuseMats = true;
}


NormBase* DarcyAdvection::getNormIntegrand (AnaSol* asol) const
{
  return new DarcyAdvectionNorm(*const_cast<DarcyAdvection*>(this));
}


bool DarcyAdvectionNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                  const TimeDomain& time, const Vec3& X) const
{
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const DarcyAdvection& problem = static_cast<const DarcyAdvection&>(myProblem);
  const double D = problem.mat->getDispersivity(X);

  // Evaluate the concentration field gradient
  Vector dCh;
  fe.grad(1).multiply(elmInt.vec.front(), dCh, true);

  pnorm[0] += D*dCh*dCh*fe.detJxW;

  return true;
}


size_t DarcyAdvectionNorm::getNoFields (int group) const
{
  return group < 1 ? 1 : (group > 1 ? 0 : 1);
}


std::string DarcyAdvectionNorm::getName (size_t, size_t, const char* prefix) const
{
  static const char* name = "concentration gradient norm";

  return prefix ? prefix + std::string(" ") + name : name;
}


