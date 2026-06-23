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
#include "FiniteElement.h"


DarcyBase::DarcyBase (unsigned short int n, int torder)
  : HasGravityBase(n), bdf(torder)
{
  ownerSim = nullptr;
  mat = nullptr;
}


Vec3 DarcyBase::getPermeability (const Vec3& X) const
{
  return mat ? mat->getPermeability(X) : Vec3();
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
    primsol.resize(1+bdf.getActualOrder());
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
