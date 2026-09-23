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

#include "SIMsolution.h"

class DarcyAdvection;
class DataExporter;
class TimeStep;


/*!
  \brief Driver class for isogeometric FE analysis of Darcy advection problems.
*/

template<class Dim>
class SIMDarcyAdvection : public Dim, public SIMsolution
{
public:
  //! \brief Default constructor.
  //! \param itg Integrand to use
  explicit SIMDarcyAdvection(DarcyAdvection& itg);

  //! \brief Destructor.
  virtual ~SIMDarcyAdvection();

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "DarcyAdvection"; }

  //! \brief Register fields for data export.
  void registerFields(DataExporter& exporter);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Initialize simulator.
  bool init();

  //! \brief Initialize time-dependent simulator.
  bool init(const TimeStep&) { return init(); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(const TimeStep& tp, bool forceNewTangent = false);

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
    return this->restoreSolution(data,this->getName());
  }

  //! \brief Initializes the material parameters for current patch.
  bool initMaterial(size_t propInd) override;

private:
  DarcyAdvection& drc; //!< Reference to the Darcy advection integrand

  std::vector<DarcyMaterial> mVec; //!< Vector of patchwise material data

  bool newTangent = true; //!< True to assemble system matrix
};

#endif
