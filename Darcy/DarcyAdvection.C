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
#include "DarcyMaterial.h"

#include "ElementSteps.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "EqualOrderOperators.h"
#include "ExprFunctions.h"
#include "Field.h"
#include "FiniteElement.h"
#include "IFEM.h"
#include "SIMbase.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "tinyxml2.h"


DarcyAdvection::DarcyAdvection (unsigned short int n, int torder)
  : DarcyBase(n,torder)
{
  this->registerVector("pressure",&pVec);

  bodyforce = nullptr;
}


DarcyAdvection::~DarcyAdvection ()
{
  delete bodyforce;
}


bool DarcyAdvection::parse (const tinyxml2::XMLElement* elem)
{
  const char* input = nullptr;
  if (!bodyforce && (input = utl::getValue(elem,"bodyforce")))
    bodyforce = new VecFuncExpr(input);
  else if ((input = utl::getValue(elem,"source")))
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
      src = new DiracSum(input,tol,nsd);
    }
    else if (type == "elementsum" && ownerSim->createFEMmodel('y'))
      src = new ElementSteps(input,*ownerSim,nsd);

    if (src)
      source.reset(src);
  }
  else
    return this->DarcyBase::parse(elem);

  return true;
}


bool DarcyAdvection::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  using WeakOps = EqualOrderOperators::Weak;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A.empty() && calcMats)
  {
    const Vec3 K = mat->getPermeability(X);
    const double D = mat->getDispersivity(X);
    const double mu = mat->getViscosity();

    // Evaluate the Darcy velocity, q = -K/mu * (grad(p) - rho*(g + bf))
    Vector dP;
    pField->gradFE(fe,dP);
    Vec3 q(dP), bf(gravity);
    if (bodyforce)
      bf += (*bodyforce)(X);
    if (!bf.isZero())
      q -= mat->getDensity(elmInt.vec.front().dot(fe.N))*bf;
    for (size_t i = 0; i < nsd; i++)
      q[i] *= -K[i]/mu;

    WeakOps::Laplacian(elMat.A.front(), fe, D, false);
    WeakOps::Advection(elMat.A.front(), fe, q, 1.0, WeakOperators::CONSERVATIVE);
  }

  if (source)
    WeakOps::Source(elMat.b.front(), fe, (*source)(X), 1);

  if (bdf.getActualOrder() > 0)
  {
    const double phi = mat->getPorosity(X);

    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++)
      c -= elmInt.vec[t].dot(fe.N) * phi * bdf[t] / time.dt;

    WeakOps::Source(elMat.b.front(), fe, c, 1);
    if (!elMat.A.empty() && calcMats)
      WeakOps::Mass(elMat.A.front(), fe, phi*bdf[0] / time.dt);
  }

  return true;
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
  s.push_back(dCh.begin(),dCh.end());

  return true;
}


void DarcyAdvection::setNamedField (const std::string& name, Field* field)
{
  if (name == "pressure")
    pField.reset(field);
  else
    delete field;
}


NormBase* DarcyAdvection::getNormIntegrand (AnaSol*) const
{
  return new DarcyAdvectionNorm(*const_cast<DarcyAdvection*>(this));
}


bool DarcyAdvectionNorm::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
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
