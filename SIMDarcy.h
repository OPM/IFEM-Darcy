// $Id$
//==============================================================================
//!
//! \file SIMDarcy.h
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Simulation driver for Isogeometric FE analysis of Darcy Flow.
//!
//==============================================================================

#ifndef _SIM_DARCY_H_
#define _SIM_DARCY_H_

#include "IFEM.h"
#include "Darcy.h"
#include "DarcySolutions.h"
#include "AnaSol.h"
#include "Functions.h"
#include "ExprFunctions.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "Property.h"
#include "DataExporter.h"


/*!
  \brief Driver class for isogeometric FE analysis of Darcy flow problems.
*/

template<class Dim> class SIMDarcy : public Dim
{
public:
  //! \brief Default constructor.
  SIMDarcy() : Dim(1), drc(Dim::dimension), solVec(&sol)
  {
    Dim::myProblem = &drc;
    aCode[0] = aCode[1] = 0;
  }

  //! \brief Destructor.
  virtual ~SIMDarcy()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"darcy"))
      return this->Dim::parse(elem);

    const char* value = nullptr;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())

      if ((value = utl::getValue(child,"permvalues")))
        drc.setPermValues(new VecFuncExpr(value));
      else if ((value = utl::getValue(child,"permeability"))) {
        std::string type;
        utl::getAttribute(child,"type",type);
        IFEM::cout <<"\tPermeability";
        drc.setPermField(utl::parseRealFunc(value,type));
        IFEM::cout << std::endl;
      }
      else if ((value = utl::getValue(child,"bodyforce")))
        drc.setBodyForce(new VecFuncExpr(value));
      else if (!strcasecmp(child->Value(),"source")) {
        std::string type;
        utl::getAttribute(child,"type",type);
        IFEM::cout <<"\tSource function: ";
        if (type == "Wavefront") {
          IFEM::cout << "Wavefront"<< std::endl;
          drc.setSource(new WavefrontSource());
        }
        else if (type == "expression" && child->FirstChild()) {
          IFEM::cout << child->FirstChild()->Value() << std::endl;
          drc.setSource(new EvalFunction(child->FirstChild()->Value()));
        }
        else
          IFEM::cout <<"(none)"<< std::endl;
      }
      else if (!strcasecmp(child->Value(),"anasol")) {
        std::string type;
        utl::getAttribute(child,"type",type);
        if (type == "Lshape") {
          Dim::mySol = new AnaSol(new LshapeDarcy(), new LshapeDarcyVelocity());
          IFEM::cout <<"\tAnalytical solution: Lshape"<< std::endl;
        }
        else if (type == "Wavefront") {
          Dim::mySol = new AnaSol(new Wavefront(), new WavefrontVelocity());
          IFEM::cout <<"\tAnalytical solution: Wavefront"<< std::endl;
        }
        else {
          Dim::mySol = new AnaSol(child);
          std::cout <<"\tAnalytical solution: expression"<< std::endl;
        }
        // Define the analytical boundary traction field
        int code = 0;
        if (utl::getAttribute(child,"code",code) && code > 0) {
          if (code > 0 && Dim::mySol->getScalarSecSol())
          {
            this->setPropertyType(code,Property::NEUMANN);
            Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
            aCode[1] = code;
          }
        }
      }
      else
        this->Dim::parse(child);

    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

    if (sit != Dim::myScalars.end())
      drc.setFlux(sit->second);
    else if (vit != Dim::myVectors.end())
      drc.setFlux(vit->second);
    else
      return false;

    return true;
  }

  //! \brief Initializes the property containers of the model.
  //! \details Use this method to clear the model before re-reading
  //! the input file in the refinement step of an adaptive simulation.
  virtual void clearProperties()
  {
    // To prevent SIMbase::clearProperties deleting the analytical solution
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
    aCode[0] = aCode[1] = 0;

    drc.setFlux((RealFunc*)nullptr);
    drc.setFlux((VecFunc*)nullptr);
    this->Dim::clearProperties();
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "DarcyFlow"; }

  //! \brief Set solution vector used.
  //! \details Used to supply an external solution vector for adaptive simulations.
  //! \param sol Pointer to vector to use.
  void setSol(const Vector* sol) { solVec = sol; }

  //! \brief Return solution vector.
  const Vector& getSolution(int=0) { return *solVec; }

  //! \brief Register fields for data export
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("p", "primary", DataExporter::SIM,DataExporter::PRIMARY);
    exporter.setFieldValue("p", this, solVec);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep&, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    // Write solution fields
    if (!this->writeGlvS(*solVec,1,nBlock,0.0,0))
      return false;

    if (!solVec->empty()) {
      Matrix tmp;
      if (!this->project(tmp, *solVec))
        return false;
      this->writeGlvV(tmp,"velocity",1,nBlock,110,Dim::nsd);
    }

    return this->writeGlvStep(1,0.0,1);
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep&)
  {
    if (!this->setMode(SIM::DYNAMIC))
      return false;

    this->initSystem(Dim::opt.solver,1,1,0,true);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);

    if (!this->assembleSystem())
      return false;

    if (!this->solveSystem(sol, Dim::msgLevel-1,"pressure    "))
      return false;

    this->printSummary();
    return true;
  }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  virtual void preprocessA()
  {
    if (!Dim::mySol) return;

    // Define analytical boundary condition fields
    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); ++p)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        if (!Dim::mySol->getScalarSol())
          p->pcode = Property::UNDEFINED;
        else if (aCode[0] == abs(p->pindx))
          p->pcode = Property::DIRICHLET_INHOM;
        else if (aCode[0] == 0)
        {
          aCode[0] = abs(p->pindx);
          Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
          p->pcode = Property::DIRICHLET_INHOM;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        if (!Dim::mySol->getScalarSecSol())
          p->pcode = Property::UNDEFINED;
        else if (aCode[1] == p->pindx)
          p->pcode = Property::NEUMANN;
        else if (aCode[1] == 0)
        {
          aCode[1] = p->pindx;
          Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
          p->pcode = Property::NEUMANN;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
  }

  //! \brief Print norm summary output.
  void printSummary()
  {
    // Evaluate solution norms
    Matrix eNorm;
    Vectors gNorm;
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (this->solutionNorms(this->getSolution(),eNorm,gNorm))
      this->printNorms(gNorm);
  }

private:
  Darcy drc;            //!< Darcy integrand.
  const Vector* solVec; //!< Pointer to solution vector.
  Vector sol;           //!< Internal solution vector.
  int aCode[2];         //!< Analytical BC code (used by destructor)
};

#endif
