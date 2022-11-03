// $Id$
//==============================================================================
//!
//! \file SIMDarcyAdvection.h
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of Darcy advection.
//!
//==============================================================================

#ifndef _SIM_DARCY_ADVECTION_H_
#define _SIM_DARCY_ADVECTION_H_

#include "DarcyMaterial.h"

#include "MatVec.h"
#include "SIMconfigure.h"
#include "SIMenums.h"
#include "SIMsolution.h"
#include "TextureProperties.h"

#include <cstddef>
#include <iosfwd>
#include <string>


class DarcyAdvection;
class DataExporter;
class TimeStep;
class TiXmlElement;
class VTF;


/*!
  \brief Driver class for isogeometric FE analysis of Darcy advection problems.
*/

template<class Dim>
class SIMDarcyAdvection : public Dim, public SIMsolution
{
public:
  //! \brief Setup properties.
  struct SetupProps {
    DarcyAdvection* itg = nullptr; //!< Pointer to integrand
  };

  //! \brief Default constructor.
  //! \param torder Order of BDF time stepping
  explicit SIMDarcyAdvection(DarcyAdvection& itg);

  //! \brief Construct from setup properties.
  //! \param torder Order of BDF time stepping
  explicit SIMDarcyAdvection(const SetupProps& props) : SIMDarcyAdvection(*props.itg) {}

  //! \brief Destructor.
  virtual ~SIMDarcyAdvection();

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "DarcyAdvection"; }

  //! \brief Register fields for data export.
  void registerFields(DataExporter& exporter);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Initialize simulator.
  bool init();

  bool init(const TimeStep&) { return init(); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(const TimeStep& tp);

  //! \brief Post-process solution.
  void postSolve(const TimeStep&) {}

  //! \brief Advance time stepping
  bool advanceStep(TimeStep&);

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override
  {
    return this->saveSolution(data,this->getName());
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    return true;
  }

  //! \brief Initializes the material parameters for current patch.
  bool initMaterial(size_t propInd) override;

  bool newTangent = true; //!< True to assemble system matrix

protected:
  DarcyAdvection& drc;           //!< Darcy integrand

private:
  std::vector<DarcyMaterial> mVec; //!< Vector of patchwise material data
  TextureProperties props; //!< Textured properties
};


//! \brief Partial specialization for configurator.
template<class Dim>
struct SolverConfigurator<SIMDarcyAdvection<Dim>> {
  //! \brief Configure a SIMDarcyAdvection instance.
  //! \param ad The SIMDarcy instance to configure
  //! \param[in] props Configuration properties
  //! \param[in] infile The input file to read
  int setup(SIMDarcyAdvection<Dim>& ad,
            const typename SIMDarcyAdvection<Dim>::SetupProps& props,
            char* infile);
};

#endif
