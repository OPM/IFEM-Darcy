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

#include "AnaSol.h"
#include "SIMDarcyTransportCorr.h"

#include "DarcyTransportCorr.h"

#include "ASMbase.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "Profiler.h"
#include "SIM1D.h"
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
SIMDarcyTransportCorr<Dim>::SIMDarcyTransportCorr (DarcyTransportCorr& itg) :
  Dim(Dim::dimension)
{
  Dim::myProblem = &itg;
  Dim::myHeading = "Darcy transport correction solver";
}


template<class Dim>
SIMDarcyTransportCorr<Dim>::
SIMDarcyTransportCorr (DarcyTransportCorr& itg,
                       const std::vector<unsigned char>& nf) : Dim(nf)
{
  Dim::myProblem = &itg;
  Dim::myHeading = "Darcy transport correction solver";
}


template<class Dim>
SIMDarcyTransportCorr<Dim>::~SIMDarcyTransportCorr ()
{
  if (vCode > 0)
    Dim::myVectors.erase(vCode);
  Dim::myProblem = nullptr;
  Dim::myInts.clear();
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"darcycorrection"))
    return this->Dim::parse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"constrain_integrated_multiplier")) {
      constrainIntegratedLag = true;
      IFEM::cout << "\tConstraining integrated multiplier";
    } else if (!strcasecmp(child->Value(), "anasol")) {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "expression")
        this->mySol = new DarcyTCorr(child);
    } else if (!Dim::myProblem->parse(child))
      this->Dim::parse(child);

  return true;
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::preprocessBeforeAsmInit (int& nnod)
{
  if (constrainIntegratedLag) {
    ++nnod;
    for (int p = 1; p <= this->getNoPatches(); ++p)
      if (int lp = this->getLocalPatchIndex(p); lp > 0)
        this->getPatch(lp)->addGlobalLagrangeMultipliers({nnod}, 1);
  }

  return true;
}


template<class Dim>
void SIMDarcyTransportCorr<Dim>::preprocessA ()
{
  Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

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
bool SIMDarcyTransportCorr<Dim>::saveModel (char* fileName, int& geoBlk, int& nBlock)
{
  return Dim::opt.format < 0 ? true : this->writeGlvG(geoBlk,fileName);
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  if (Dim::opt.format < 0 || (tp.step % Dim::opt.saveInc) > 0)
    return true;

  // Write solution fields
  if (!this->writeGlvS(qSol,1,nBlock,0.0,"Darcy velocity"))
    return false;

  // Write element norms
  if (Dim::opt.saveNorms)
    if (!this->writeGlvN(eNorm,1,nBlock))
      return false;

  return this->writeGlvStep(1, 0.0, 1);
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::init ()
{
  this->initSystem(Dim::opt.solver,1,1);
  this->setQuadratureRule(Dim::opt.nGauss[0],true);
  return true;
}


template<class Dim>
bool SIMDarcyTransportCorr<Dim>::solveStep (const TimeStep&)
{
  if (!this->setMode(SIM::DYNAMIC))
    return false;
  if (!this->initDirichlet())
    return false;
  if (!this->assembleSystem())
    return false;

  return this->solveSystem(qSol,Dim::msgLevel-1,"darcy velocity");
}


template<class Dim>
void SIMDarcyTransportCorr<Dim>::
printSolutionSummary (const Vector&, int, const char*, std::streamsize outPrec)
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
        str <<"\n               Max "<< char('X'+d)
            <<"-Darcy velocity : "<< dMax[d] <<" node "<< iMax[d];
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
      str <<"\nL2(q - q^) / L2(q^) = "<< diff
          <<"%\nL2(div q^) / L2(q^) = "<< divQ
          <<"%\nL2(res q^) / L2(c*q^) = "<< resQ
          <<"%\nL2(div q) / L2(q) = "<< divq
          <<"% "<< gNorm.front()[6] <<" "<< gNorm.front()[5];
  }

  IFEM::cout << str.str() << std::endl;
}


template<class Dim>
int SolverConfigurator<SIMDarcyTransportCorr<Dim>>::
setup (SIMDarcyTransportCorr<Dim>& darcy,
       const typename SIMDarcyTransportCorr<Dim>::SetupProps& props,
       char* infile)
{
  utl::profiler->start("Model input");

  if (!darcy.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (!darcy.preprocess())
    return 2;

  darcy.init();

  return 0;
}


//! \brief Instantiation macro.
#define INSTANTIATE(T) \
  template class SIMDarcyTransportCorr<T>; \
  template struct SolverConfigurator<SIMDarcyTransportCorr<T>>;

INSTANTIATE(SIM1D)
INSTANTIATE(SIM2D)
INSTANTIATE(SIM3D)
