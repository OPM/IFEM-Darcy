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

#ifndef SIMDARCY_H_
#define SIMDARCY_H_

#include "Darcy.h"
#include "ASMbase.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "Property.h"
#include "DataExporter.h"
#include "AnaSol.h"
#include "DarcySolutions.h"


template<class Dim> class SIMDarcy : public Dim
{
  public:
    SIMDarcy() : Dim(1), DRC(Dim::dimension)
    {
      Dim::myProblem = &DRC;
      aCode[0] = aCode[1] = 0;
    }

    virtual ~SIMDarcy()
    {
      Dim::myProblem = NULL;
      // To prevent the SIMbase destructor try to delete already deleted functions
      if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
      if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
    }

    virtual bool parse(const TiXmlElement* elem)
    {
      if (strcasecmp(elem->Value(),"darcy"))
        return this->Dim::parse(elem);

      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement()) {
        if (!strcasecmp(child->Value(),"permvalues")) {
          const char* value = utl::getValue(child, "permvalues");
          if (value) {
            std::stringstream str;
            str << value;
            DRC.setPermValues(new VecFuncExpr(value));
          }
        }
        if (!strcasecmp(child->Value(),"permeability")){
          std::string type;
          utl::getAttribute(child, "type", type);
          RealFunc* permeability = utl::parseRealFunc(utl::getValue(child, "permeability"), type);
          std::cout << "permeability = " << permeability << std::endl;
          DRC.setPermField(permeability);
        }
        if (!strcasecmp(child->Value(),"bodyforce")) {
          const char* value = utl::getValue(child, "bodyforce");
          if (value) {
            std::stringstream str;
            str << value;
            DRC.setBodyForce(new VecFuncExpr(value));
          }
        }
        if (!strcasecmp(child->Value(),"source")) {
          std::string type;
          utl::getAttribute(child,"type",type);
          if (type == "Wavefront") {
            std::cout << "Source function: Wavefront" << std::endl;
            DRC.setSource(new WavefrontSource());
          }
          else if (type == "expression" && child->FirstChild()) {
            std::cout << "Source function: " << child->FirstChild()->Value() << std::endl;
            DRC.setSource(new EvalFunction(child->FirstChild()->Value()));
          }
        }
        else if (!strcasecmp(child->Value(),"anasol")) {

          std::string type;
          utl::getAttribute(child,"type",type);
          if (type == "Lshape") {
            this->mySol = new AnaSol(new LshapeDarcy(), new LshapeDarcyVelocity());
            std::cout << "Anasol: Lshape" << std::endl;
          }
          else if (type == "Wavefront") {
            this->mySol = new AnaSol(new Wavefront(), new WavefrontVelocity());
            std::cout << "Anasol: Wavefront" << std::endl;
          } else {
            this->mySol = new AnaSol(child);
            std::cout << "Anasol: expression" << std::endl;
          }
          // Define the analytical boundary traction field
          int code = 0;
          if (code == 0 && utl::getAttribute(child,"code",code)) {
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
      }
      return true;
    }

    virtual bool initNeumann(size_t propInd)
    {
      typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
      typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

      if (sit != Dim::myScalars.end())
        DRC.setFlux(sit->second);
      else if (vit != Dim::myVectors.end())
        DRC.setFlux(vit->second);
      else
        return false;

      return true;
    }

    virtual void clearProperties()
    {
      // To prevent SIMbase::clearProperties deleting the analytical solution
      if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
      if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
      aCode[0] = aCode[1] = 0;

      DRC.setFlux((RealFunc*)NULL);
      DRC.setFlux((VecFunc*)NULL);
      this->Dim::clearProperties();
    }

    virtual std::string getName() const
    {
      return "DarcyFlow";
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
      {
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
    }

  private:
    Darcy DRC;
    Vectors pressure;
    int aCode[2];   //!< Analytical BC code (used by destructor)
};

#endif /* SIMDARCY_H_ */
