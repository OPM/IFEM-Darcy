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

#include "Darcy.h"
#include "SIMsolution.h"


class DataExporter;
class TimeStep;


/*!
  \brief Driver class for isogeometric FE analysis of Darcy flow problems.
*/

template<class Dim> class SIMDarcy : public Dim, public SIMsolution
{
public:
  //! \brief Default constructor.
  SIMDarcy(int torder = 0);

  //! \brief Destructor.
  virtual ~SIMDarcy();

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override;

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

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
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Initialize simulator.
  void init();

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp);

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
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  //!
  //! This overloaded version also computes the reaction forces along a given
  //! boundary. This requires an additional assembly loop calculating the
  //! internal forces only, since we only are doing a linear solve here.
  bool solveSystem(Vector& solution, int printSol, double* rCond,
                   const char* compName, bool newLHS, size_t idxRHS) override;

  //! \brief Returns current reaction force vector.
  const Vector* getReactionForces() const override
  {
    return myReact.empty() ? nullptr : &myReact;
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data)
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
  //\param tp Time stepping parameters
  void printFinalNorms(const TimeStep& tp);

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  void preprocessA() override;

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocessB() override;

private:
  Darcy drc;            //!< Darcy integrand
  const Vector* solVec; //!< Pointer to solution vector
  Vector myReact;       //!< Nodal reaction forces
  int aCode[2];         //!< Analytical BC code (used by destructor)
  Matrix eNorm;         //!< Element wise norms
  Vectors proj;         //!< Projected solution vectors
};

#endif
