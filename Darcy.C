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
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "AnaSol.h"
#include "WeakOperators.h"


Darcy::Darcy(unsigned short int n) :
  flux(NULL), source(NULL), vflux(NULL), nsd(n), rhow(1.0), gacc(9.81)
{
  primsol.resize(1);
  source = 0;
  bodyforce = 0;
  permeability = 0;
  permvalues = 0;
}


double Darcy::getPotential (const Vec3& X) const
{
  return source ? (*source)(X) : 0.0;
}


double Darcy::getFlux(const Vec3& X, const Vec3& normal) const
{
  if (flux)
    return (*flux)(X);
  else if (vflux)
    return (*vflux)(X)*normal;
  else
    return 0.0;
}


LocalIntegral* Darcy::getLocalIntegral(size_t nen, size_t,
                                            bool neumann) const
{
  ElmMats* result = new ElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(neumann ? 0 : 1, 1);


  if (!neumann)
    result->A[0].resize(nen, nen, true);
  result->b[0].resize(nen, true);
  result->redim(nen);

  return result;
}


bool Darcy::evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                    const Vec3& X) const
{

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // integration for permeability matrix and the rhs vector
  if (!elMat.A.empty())
  {
    // Evaluate the hydraulic conductivity matrix at this point
    Matrix K;
    this->formKmatrix(K,X);
    WeakOperators::LaplacianCoeff(elMat.A.front(), K, fe, 1.0/(rhow*gacc));
  }

  Vec3 b = bodyforce ? (*bodyforce)(X) : Vec3();
  // Read components of permeability matrix
  Vec3 perm = permvalues ? (*permvalues)(X) : Vec3();
  // Read point-wise permeabilities from a random field
  double kf = permeability ? (*permeability)(X) : 0.0;

  Vec3 eperm = b;
  if (permvalues)
    for (size_t i=0;i<nsd;++i)
      eperm[i] *= perm[i]/gacc;

  if (permeability)
    eperm *= kf/gacc;

  // Integrate rhs contribution from body force
  WeakOperators::Divergence(elMat.b[0], fe, eperm);

  if (source)
    WeakOperators::Source(elMat.b.front(), fe, (*source)(X));

  return true;
}


bool Darcy::evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                    const Vec3& X, const Vec3& normal) const
{
  if (!flux && !vflux)
  {
    std::cerr << " *** Darcy::evalBou: No fluxes." << std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty())
  {
    std::cerr << "*** Darcy::evalBou: No load vector." << std::endl;
    return false;
  }

  double qw = -this->getFlux(X,normal);

  WeakOperators::Source(elMat.b.front(), fe, qw/rhow);

  return true;
}


bool Darcy::formKmatrix (Matrix& K, const Vec3& X, bool inverse) const
{
  Vec3 perm = permvalues ? (*permvalues)(X) : Vec3();

  double kf = permeability ? (*permeability)(X) : 0.0;

  K.resize(nsd,nsd);
  for (int i = 1; i <= nsd; i++)
  {
    if (permvalues)
      K(i,i) = perm[i-1];
    if (permeability)
      K(i,i) = kf;
  }

  if (inverse)
    K.inverse();

  return true;
}


bool Darcy::evalSol (Vector& v, const FiniteElement& fe,
                      const Vec3& X, const std::vector<int>& MNPC) const
{
  if (primsol.empty() || primsol.front().empty())
  {
    std::cerr << "*** Darcy::evalSol: No primary solution." << std::endl;
    return false;
  }

  Vector eV;
  int ierr = utl::gather(MNPC,1,primsol.front(),eV);
  if (ierr > 0)
  {
    std::cerr << "*** Darcy::evalSol: Detected " << ierr
              << " node numbers out of range." <<  std::endl;
    return false;
  }

  return this->evalSol(v,eV,fe.dNdX,X);
}


bool Darcy::evalSol (Vector& v, const Vector& eV,
                     const Matrix& dNdX, const Vec3& X) const
{
  if (eV.size() != dNdX.rows())
  {
    std::cerr << "*** Darcy:evalSol: Invalid solution vector."
              << "\n  size(eV) = " << eV.size() << " size(dNdX) = "
              << dNdX.rows() << "," << dNdX.cols() << std::endl;
    return false;
  }

  Matrix K,KB;
  this->formKmatrix(K,X);

  Vec3 b = bodyforce ? (*bodyforce)(X) : Vec3();

  std::vector<double> temp(v.size());

  // Evaluate the velocity
  dNdX.multiply(eV, temp, true);
  for (size_t i=0;i<nsd;++i)
    temp[i] -= rhow*b[i];

  K.multiply(temp,v);
  v *= -1.0/(rhow*gacc);

  return true;
}


bool Darcy::evalSol(Vector& s, const VecFunc& asol,
                                 const Vec3& X) const
{
  s = Vector(asol(X).ptr(),nsd);
  return true;
}

const char* Darcy::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "pressure";

  static std::string name;
  name = prefix + std::string(" pressure");

  return name.c_str();
}


const char* Darcy::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return 0;

  static const char* s[3] = {"v_x","v_y","v_z"};
  if(!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


NormBase* Darcy::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new DarcyNorm(*const_cast<Darcy*>(this),
                           asol->getScalarSecSol());
  else
    return new DarcyNorm(*const_cast<Darcy*>(this));
}


Vec3 Darcy::getPermeability(const Vec3& X) const
{
  Vec3 result;
  if (permvalues)
    result = (*permvalues)(X);
  else if (permeability)
    result = (*permeability)(X);

  return result;
}


Vec3 Darcy::getBodyForce(const Vec3& X) const
{
  Vec3 result;
  if (bodyforce)
    result = (*bodyforce)(X);

  return result;
}


DarcyNorm::DarcyNorm (Darcy& p, VecFunc* a) : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


bool DarcyNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X) const
{
  Darcy& problem = static_cast<Darcy&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Kinv;
  problem.formKmatrix(Kinv,X,true);
  double rg = (problem.rhow)*(problem.gacc);

  // Evaluate the finite element pressure field
  Vector sigmah, sigma, error;

  //fe.dNdX.multiply(pnorm.vec.front(),sigmah,true);
  if (!problem.evalSol(sigmah,pnorm.vec.front(),fe.dNdX,X))
    return false;

  size_t ip = 0;
  // Integrate the energy norm a(u^h,u^h)
  pnorm[ip++] += sigmah.dot(Kinv*sigmah)*rg*fe.detJxW;
  // Evaluate the pressure field
  double p = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (h,u^h)
  pnorm[ip++] += problem.getPotential(X)*p*fe.detJxW;

  if (anasol)
  {
    // Evaluate the analytical velocity
    sigma.fill((*anasol)(X).ptr(),nrcmp);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Kinv*sigma)*rg*fe.detJxW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Kinv*error)*rg*fe.detJxW;
  }

  size_t i,j;

  for (i = 0; i < pnorm.psol.size(); i++)
  {
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected heat flux field
      Vector sigmar(nrcmp);

      for (j = 0; j < nrcmp; j++)
        sigmar[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Kinv*sigmar)*rg*fe.detJxW;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Kinv*error)*rg*fe.detJxW;

      if (anasol)
      {
        // Integrate the error in the projected solution a(u-u^r,u-u^r)
        error = sigma - sigmar;
        pnorm[ip++] += error.dot(Kinv*error)*rg*fe.detJxW;
        ip++; // Make room for the local effectivity index here
      }
    }
  }

  return true;
}


bool DarcyNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X, const Vec3& normal) const
{
  Darcy& problem = static_cast<Darcy&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface heat flux
  double h = problem.getFlux(X,normal);
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);

  // Integrate the external energy (h,u^h)
  pnorm[1] += h*u*fe.detJxW;
  return true;
}


bool DarcyNorm::finalizeElement (LocalIntegral& elmInt,
                                 const TimeDomain&, size_t)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 7; ip < pnorm.size(); ip += 4)
    pnorm[ip] = sqrt(pnorm[ip-2] / pnorm[3]);

  return true;
}


void DarcyNorm::addBoundaryTerms (Vectors& gNorm, double energy) const
{
  gNorm.front()[1] += energy;
}


size_t DarcyNorm::getNoFields (int group) const
{
  size_t nf = 1;
  if (group < 1)
    for (size_t i = 0; i < prjsol.size(); i++)
      nf += prjsol.empty() ? 0 : 1;
  else
    nf = anasol ? 4 : 2;

  return nf;
}


const char* DarcyNorm::getName (size_t i, size_t j, const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[8] = {
    "a(u^h,u^h)^0.5",
    "(h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index"
  };

  size_t k = i > 1 ? j+3 : j-1;

  if (!prefix)
    return s[k];

  static std::string name;
  name = prefix + std::string(" ");
  name += s[k];

  return name.c_str();
}


bool DarcyNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}
