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

#include "DarcyEnums.h"
#include "DarcyMaterial.h"

#include "MatVec.h"
#include "SIMconfigure.h"
#include "SIMenums.h"
#include "SIMMultiPatchModelGen.h"
#include "SIMsolution.h"

#include <cstddef>
#include <iosfwd>
#include <string>


class Darcy;
class DataExporter;
class TimeStep;
class VTF;


/*!
  \brief Driver class for isogeometric FE analysis of Darcy flow problems.
*/

template<class Dim>
class SIMDarcy : public SIMMultiPatchModelGen<Dim>, public SIMsolution
{
public:
  //! \brief Setup properties.
  struct SetupProps {
    Darcy* itg = nullptr; //!< Pointer to integrand
  };

  //! \brief Default constructor.
  //! \param torder Order of BDF time stepping
  //! \param nf Number of primary fields
  explicit SIMDarcy(Darcy& itg, unsigned char nf = 1);

  //! \brief Default constructor.
  //! \param torder Order of BDF time stepping
  //! \param nf Number of primary fields
  SIMDarcy(Darcy& itg, const std::vector<unsigned char>& nf);

  //! \brief Construct from setup properties.
  //! \param torder Order of BDF time stepping
  explicit SIMDarcy(const SetupProps& props) : SIMDarcy(*props.itg) {}

  //! \brief Destructor.
  virtual ~SIMDarcy();

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

  //! \brief Initializes the material parameters for current patch.
  bool initMaterial(size_t propInd) override;

  //! \brief Initializes the property containers of the model.
  //! \details Use this method to clear the model before re-reading
  //! the input file in the refinement step of an adaptive simulation.
  void clearProperties() override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "DarcyFlow"; }

  //! \brief Set solution vector used.
  //! \details Used to supply an external solution vector for adaptive simulations.
  //! \param sol Pointer to vector to use.
  void setSol(const Vector* sol) { solVec = sol; }

  //! \brief Return solution vector.
  const Vector& getSolution(int = 0) const override { return *solVec; }

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
  //! \param[in] newData If true we have new pressure data
  bool saveStep(const TimeStep& tp, int& nBlock, bool newData = true);

  //! \brief Initialize simulator.
  bool init();

  bool init(const TimeStep&) { return init(); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(const TimeStep& tp);

  //! \brief Post-process solution.
  void postSolve(const TimeStep&) {}

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  //!
  //! \details Since this solver is linear, this is just a normal solve.
  SIM::ConvStatus solveIteration(TimeStep& tp);

  //! \brief Advance time stepping
  bool advanceStep(TimeStep&);

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] printSol Print solution only if size is less than this value
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solution, int printSol, const char*,
                            std::streamsize outPrec) override;

  using Dim::solveSystem;
  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  //!
  //! \details This overloaded version also computes the reaction forces along
  //! a given boundary. This requires an additional assembly loop calculating
  //! the internal forces only, since we only are doing a linear solve here.
  bool solveSystem(Vector& solution, int printSol, double* rCond,
                   const char* compName, size_t idxRHS) override;

  //! \brief Returns current reaction force container.
  const RealArray* getReactionForces() const override
  {
    return myReact.empty() ? nullptr : &myReact;
  }

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

  //! \brief Print final solution norms to terminal.
  //! \param tp Time stepping parameters
  void printFinalNorms(const TimeStep& tp);

  //! \brief Print solution solution norms to terminal.
  void printSolNorms(const Vector& gNorm, size_t w) const;

  //! \brief Print norms to screen during adaptive simulations.
  void printExactNorms (const Vector& gNorm, size_t w = 36) const;

  //! \brief Print norms to screen during adaptive simulations.
  void printNorms (const Vectors& gNorm, size_t w = 36) const override;

  //! \brief Prints a norm group to the log stream.
  void printNormGroup(const Vector& rNorm, const Vector& fNorm,
                      const std::string& name) const override;

  using Dim::savePoints;
  //! \brief Saves point results to output file for a given time step.
  //! \param[in] time Load/time step parameter
  //! \param[in] iStep Load/time step counter
  bool savePoints(double time, int iStep) const
  {
    return this->savePoints(*solVec, time, iStep);
  }

  //! \brief Set norm to adapt based on.
  //! \param norm Norm to use
  void setAdaptiveNorm(DCY::AdaptationNorm norm)
  {
    adNorm = norm;
  }

  //! \brief Returns norm to adapt based on.
  DCY::AdaptationNorm getAdaptiveNorm() const
  {
    return adNorm;
  }

  //! \brief Returns the reference norm to base mesh adaptation upon.
  //! \param[in] gNorm The calculated global norms
  //! \param[in] adaptor Which norm group to base adaptation on
  double getReferenceNorm(const Vectors& gNorm, size_t adaptor) const override;

  //! \brief Returns the global effectivity index.
  //! \param[in] gNorm Global norm values
  //! \param[in] idx 0-based norm group index
  //! \param[in] inorm 1-based norm index within the specified norm group
  double getEffectivityIndex(const Vectors& gNorm,
                             size_t adaptor,
                             size_t inorm) const override;

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  void preprocessA() override;

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocessB() override;

  //! \brief Perform solution projection.
  bool doProjection();

  Darcy& drc;           //!< Darcy integrand

private:
  DCY::AdaptationNorm adNorm = DCY::NO_ADAP; //!< Norm to adapt based on
  const Vector* solVec; //!< Pointer to solution vector
  RealArray myReact;    //!< Nodal reaction forces
  int aCode[2];         //!< Analytical BC code (used by destructor)
  Matrix eNorm;         //!< Element wise norms
  Vectors proj;         //!< Projected solution vectors

  std::vector<DarcyMaterial> mVec; //!< Vector of patchwise material data

  int maxCycle = -1;      //!< Max number of sub-iterations
  double cycleTol = 1e-6; //! < Convergence tolerance in sub-iterations

  bool newTangent = true; //!< True to assemble element matrices
};


//! \brief Partial specialization for configurator.
template<class Dim>
struct SolverConfigurator<SIMDarcy<Dim>> {
  //! \brief Configure a SIMDarcy instance.
  //! \param ad The SIMDarcy instance to configure
  //! \param[in] props Configuration properties
  //! \param[in] infile The input file to read
  int setup(SIMDarcy<Dim>& ad,
            const typename SIMDarcy<Dim>::SetupProps& props,
            char* infile);
};

#endif
