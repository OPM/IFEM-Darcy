// $Id$
//==============================================================================
//!
//! \file SIMDarcyTransportCorr.C
//!
//! \date Jan 20 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Darcy transport correction problem.
//!
//==============================================================================

#include "SIMDarcyTransportCorr.h"

#include "AlgEqSystem.h"
#include "AnaSol.h"
#include "ASMbase.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "IntegrandBase.h"
#include "Profiler.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMenums.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "tinyxml2.h"


namespace
{
  //! \brief Analytical solution class for DarcyTransportCorr
  class DarcyTCorr : public AnaSol
  {
  public:
    //! \brief Constructor initializing expression functions from XML tags.
    explicit DarcyTCorr(const tinyxml2::XMLElement* xml) : AnaSol(xml,false) {}
    //! \brief Override parent class method to avoid the secondary solution.
    void setupSecondarySolutions() override {}
  };
}


template<class Dim>
SIMDarcyTransportCorr<Dim>::SIMDarcyTransportCorr (IntegrandBase& itg,
                                                   const CharVec& nf,
                                                   bool AL) : Dim(nf)
{
  Dim::myProblem = &itg;
  Dim::myHeading = "Darcy transport correction solver";
  Dim::lagMTOK = true; // can do multithreading with global Lagrange multipliers
  useAL = AL;
  vCode = 0;
  nSclr = AL ? 3 : 0;
  maxit = 10;
  eps[0] = eps[1] = eps[2] = 1.0e99;

  if (useAL)
    Dim::msgLevel = 1;
}


template<class Dim>
SIMDarcyTransportCorr<Dim>::~SIMDarcyTransportCorr ()
{
  Dim::myVectors.erase(vCode);
  Dim::myProblem = nullptr;
  Dim::myInts.clear();
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"darcycorrection"))
    return this->Dim::parse(elem);

  utl::getAttribute(elem,"maxit",maxit);
  utl::getAttribute(elem,"eps0",eps[0]);
  utl::getAttribute(elem,"eps1",eps[1]);
  utl::getAttribute(elem,"eps2",eps[2]);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"anasol"))
    {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "expression")
        Dim::mySol = new DarcyTCorr(child);
    }
    else if (!Dim::myProblem->parse(child))
      this->Dim::parse(child);

  return true;
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::preprocessBeforeAsmInit (int& nnod)
{
  if (!this->mixedProblem())
    useAL = false; // Use penalty approach
  else if (!useAL)
  {
    // Use classical Lagrange multiplier approach with global constraints
    ++nnod;
    for (int p = 1; p <= this->getNoPatches(); ++p)
      if (int lp = this->getLocalPatchIndex(p); lp > 0)
        this->getPatch(lp)->addGlobalLagrangeMultipliers({nnod},Dim::nf[1]);
  }

  return true;
}


template<class Dim>
void SIMDarcyTransportCorr<Dim>::preprocessA ()
{
  for (Property& p : Dim::myProps)
    if (p.pcode == Property::DIRICHLET_ANASOL) {
      if (vCode == abs(p.pindx))
        p.pcode = Property::DIRICHLET_INHOM;
      else if (vCode == 0) {
        vCode = abs(p.pindx);
        Dim::myVectors[vCode] = Dim::mySol->getVectorSol();
        p.pcode = Property::DIRICHLET_INHOM;
      } else
        p.pcode = Property::UNDEFINED;
    }
}


template<class Dim>
void SIMDarcyTransportCorr<Dim>::registerFields (DataExporter& exporter)
{
  int results = DataExporter::PRIMARY;

  if (!Dim::opt.pSolOnly)
    results |= DataExporter::SECONDARY;

  if (Dim::opt.saveNorms)
    results |= DataExporter::NORMS;

  exporter.registerField("q", "primary", DataExporter::SIM, results);
  exporter.setFieldValue("q", this, &qSol, nullptr,
                         (results & DataExporter::NORMS) ? &eNorm : nullptr);
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::saveModel (char* fileName,
                                            int& geoBlk, int& nBlock)
{
  return Dim::opt.format < 0 ? true : this->writeGlvG(geoBlk,fileName);
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  if (Dim::opt.format < 0 || (tp.step % Dim::opt.saveInc) > 0)
    return true;

  const int iType = tp.multiSteps() ? 0 : 1;
  const int iStep = iType == 0 ? tp.step/Dim::opt.saveInc : 1;

  // Write solution fields
  if (!this->writeGlvS(qSol,iStep,nBlock,tp.time.t,"Darcy velocity"))
    return false;

  // Write element norms
  if (Dim::opt.saveNorms)
    if (!this->writeGlvN(eNorm,iStep,nBlock))
      return false;

  return this->writeGlvStep(iStep, iType == 1 ? iStep : tp.time.t, iType);
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::solveStep (TimeStep& tp)
{
  if (Dim::msgLevel >= 0 && tp.multiSteps())
    IFEM::cout <<"\n  step = "<< tp.step
               <<"  time = "<< tp.time.t << std::endl;

  if (!this->setMode(SIM::DYNAMIC))
    return false;

  if (!this->initDirichlet())
    return false;

  if (useAL)
    return this->solveALStep(tp);

  if (!this->assembleSystem())
    return false;

  return this->solveSystem(qSol,Dim::msgLevel,"velocity");
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::solveALStep (TimeStep& tp)
{
  bool newLHS = true;
  IFEM::cout << std::endl;

  for (tp.iter = 0; tp.iter < maxit; tp.iter++)
    if (tp.iter > 0 && this->checkConvergence(tp))
      break;
    else if (!this->assembleSystem(tp.time,{qSol},newLHS))
      return false;
    else if (!this->solveSystem(qSol))
      return false;
    else
      newLHS = false;

  this->printSolutionSummary(qSol,Dim::msgLevel,"velocity");
  return true;
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::checkConvergence (TimeStep& tp)
{
  std::vector<double> rs(nSclr);
  for (int i = 0; i < nSclr; i++)
    rs[i] = sqrt(Dim::myEqSys->getScalar(i));

  std::stringstream str;
  if (Dim::adm.getProcId() == 0)
  {
    str <<"  iter="<< tp.iter-1;
    for (size_t i = 1; i <= rs.size(); i++)
      str <<" r"<< i <<"="<< rs[i-1];
  }
  IFEM::cout << str.str() << std::endl;

  for (int i = 0; i < nSclr; i++)
    if (rs[i] > eps[i])
      return false;

  return true;
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::advanceStep (TimeStep&)
{
  return true;
}


template<class Dim>
void SIMDarcyTransportCorr<Dim>::printSolutionSummary (const Vector&, int,
                                                       const char*,
                                                       std::streamsize outPrec)
{
  const size_t nsd = this->getNoSpaceDim();

  size_t iMax[3];
  double dMax[3];
  double dNorm = this->solutionNorms(qSol,dMax,iMax,nsd);

  std::stringstream str;
  if (Dim::adm.getProcId() == 0)
  {
    if (outPrec > 0) str.precision(outPrec);

    str <<"  Primary solution summary: L2-norm : "<< utl::trunc(dNorm);

    char D = 'X';
    for (size_t d = 0; d < nsd; d++, D++)
      if (utl::trunc(dMax[d]) != 0.0)
        str <<"\n"<< std::string(20,' ') <<" Max "<< char('X'+d)
            <<"-velocity : "<< dMax[d] <<" node "<< iMax[d];
  }

  this->setQuadratureRule(Dim::opt.nGauss[1]);
  this->setMode(SIM::RECOVERY);

  // Evaluate solution norms
  Vectors gNorm;
  if (this->solutionNorms(qSol,eNorm,gNorm))
  {
    double diff = 100.0*gNorm.front()[2]/gNorm.front()[0];
    double divQ = 100.0*gNorm.front()[3]/gNorm.front()[0];
    double resQ = 100.0*gNorm.front()[4]/gNorm.front()[1];
    double divq = 100.0*gNorm.front()[6]/gNorm.front()[5];

    if (Dim::adm.getProcId() == 0)
      str <<"\n  L2(q - q^) / L2(q^)   = "<< diff <<"% of "<< gNorm.front()[0]
          <<"\n  L2(div q^) / L2(q^)   = "<< divQ <<"% of "<< gNorm.front()[0]
          <<"\n  L2(res q^) / L2(c*q^) = "<< resQ <<"% of "<< gNorm.front()[1]
          <<"\n  L2(div q ) / L2(q )   = "<< divq <<"% of "<< gNorm.front()[5];
  }

  IFEM::cout << str.str() << std::endl;
}

// Instantiations for different dimensions.

template class SIMDarcyTransportCorr<SIM2D>;
template class SIMDarcyTransportCorr<SIM3D>;
