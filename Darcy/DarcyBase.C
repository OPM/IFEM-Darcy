// $Id$
//==============================================================================
//!
//! \file DarcyBase.C
//!
//! \date Jun 23 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for the Darcy flow integrands.
//!
//==============================================================================

#include "DarcyBase.h"
#include "DarcyMaterial.h"

#include "ElmMats.h"
#include "ExprFunctions.h"
#include "FiniteElement.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml2.h"
#include <cstring>


DarcyBase::DarcyBase (unsigned short int n) : HasGravityBase(n)
{
  ownerSim  = nullptr;
  mat       = nullptr;
  bodyforce = nullptr;
}


DarcyBase::~DarcyBase ()
{
  delete bodyforce;
}


bool DarcyBase::parse (const tinyxml2::XMLElement* elem)
{
  const char* input = bodyforce ? nullptr : utl::getValue(elem,"bodyforce");
  if (!input) return this->HasGravityBase::parse(elem);

  if (strstr(input,"|"))
  {
    ExprEval::numError = 0;
    VecFunc* bf = new VecFuncExpr(input);
    if (ExprEval::numError > 0)
      delete bf;
    else
      bodyforce = bf;
  }
  else
  {
    Vec3 bf;
    char* ff = strdup(input);
    char* ss = strtok(ff," \t\r\n");
    for (int i = 0; i < 3 && ss; i++, ss = strtok(nullptr," \t\r\n"))
      bf[i] = atof(ss);
    free(ff);
    bodyforce = new ConstVecFunc(bf);
  }
  return true;
}


Vec3 DarcyBase::getBodyForce (const Vec3& X) const
{
  return bodyforce ? (*bodyforce)(X) : Vec3();
}


bool DarcyBase::getInvPermeability (const Vec3& X, Matrix& Kinv) const
{
  if (!mat) return false;

  if (mat->getPermeability(&Kinv))
  {
    if (Kinv.inverse() > 0.0)
      Kinv *= mat->getViscosity();
    else
      return false;
  }
  else
  {
    Kinv.diag(mat->getPermeability(X).vec(nsd));
    for (size_t i = 1; i <= Kinv.cols(); i++)
      if (double& k = Kinv(i,i); k != 0.0)
        k = mat->getViscosity()/k;
      else
        return false;
  }

  return true;
}


double DarcyBase::getDispersivity (const Vec3& X) const
{
  return mat ? mat->getDispersivity(X) : 0.0;
}


double DarcyBase::getViscosity () const
{
  return mat ? mat->getViscosity() : 0.0;
}


void DarcyBase::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  if (mode >= SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.resize(1+this->getOrder());
}


LocalIntegral* DarcyBase::getLocalIntegral (size_t nen, size_t,
                                            bool neumann) const
{
  ElmMats* result = new ElmMats(!neumann);

  result->rhsOnly = neumann || !calcMats || m_mode >= SIM::RHS_ONLY;
  result->resize(neumann ? 0 : 1, 1);
  result->redim(nen);

  return result;
}


void DarcyBase::initLHSbuffers (size_t nEl)
{
  if (calcMats < 0)
    return;

  if (nEl > 1)
    myKmats.resize(nEl);

  if (nEl > 0)
    calcMats = true;
  else if (!myKmats.empty())
    calcMats = false;
}


bool DarcyBase::initElement (const std::vector<int>& MNPC,
                             const FiniteElement& fe,
                             const Vec3& XC,
                             size_t nPt, LocalIntegral& elmInt)
{
  if (fe.iel > 0 && !calcMats)
  {
    size_t iel = fe.iel - 1;
    ElmMats* A = dynamic_cast<ElmMats*>(&elmInt);
    if (A && iel < myKmats.size() && !A->A.empty())
      A->A.front() = myKmats[iel];
  }

  return this->HasGravityBase::initElement(MNPC,fe,XC,nPt,elmInt);
}


bool DarcyBase::finalizeElement (LocalIntegral& elmInt,
                                 const FiniteElement& fe,
                                 const TimeDomain& time, size_t iGP)
{
  if (fe.iel > 0 && calcMats)
  {
    size_t iel = fe.iel - 1;
    ElmMats* A = dynamic_cast<ElmMats*>(&elmInt);
    if (A && iel < myKmats.size())
      myKmats[iel] = A->getNewtonMatrix();
  }

  return this->HasGravityBase::finalizeElement(elmInt,fe,time,iGP);
}
