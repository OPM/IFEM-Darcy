// $Id$
//==============================================================================
//!
//! \file SIMDarcyTransportCorr.h
//!
//! \date Aug 22 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Darcy transport correction problem.
//!
//==============================================================================

#ifndef _SIM_DARCY_TRANSPORT_CORR_H_
#define _SIM_DARCY_TRANSPORT_CORR_H_

#include "MatVec.h"

class DataExporter;
class IntegrandBase;
class TimeStep;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Driver class for Darcy transport correction problems.
*/

template<class Dim>
class SIMDarcyTransportCorr : public Dim
{
  using CharVec = std::vector<unsigned char>; //!< Convenience type alias

public:
  //! \brief Constructor.
  //! \param itg The integrand to use
  //! \param[in] nf Number of fields on each basis
  //! \param[in] AL If \e true, use the Augmented Lagrange formulation
  SIMDarcyTransportCorr(IntegrandBase& itg, const CharVec& nf, bool AL = false);
  //! \brief Destructor.
  ~SIMDarcyTransportCorr() override;

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "DarcyTransportCorr"; }

  //! \brief Returns the total number of scalar quantities to integrate.
  int getNoScalars() const { return nSclr; }

  //! \brief Register fields for data export.
  void registerFields(DataExporter& exporter);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running VTF block counter (updated)
  bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Computes the solution for the current time step.
  //! \param tp Time stepping parameters
  bool solveStep(TimeStep& tp);

  //! \brief Advances the time stepping.
  bool advanceStep(TimeStep&);

  //! \brief Prints out a summary of the calculated solution.
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector&, int, const char* compName,
                            std::streamsize outPrec = 0) override;

protected:
  //! \brief Adds two global multipliers for the constraints.
  //! \param nnod Number of nodes in the model (updated)
  bool preprocessBeforeAsmInit(int& nnod) override;

  //! \brief Performs some pre-processing tasks on the FE model.
  void preprocessA() override;

  //! \brief Solves the current time step by Augmented Lagrange.
  //! \param tp Time stepping parameters
  bool solveALStep(TimeStep& tp);

  //! \brief Checks if the Augmented Lagrange iterations have converged.
  //! \param tp Time stepping parameters
  bool checkConvergence(TimeStep& tp);

private:
  Vector qSol;  //!< Solution vector
  Matrix eNorm; //!< Element norms
  bool   useAL; //!< If \e true, use the Augmented Lagrange formulation
  int    vCode; //!< Index for analytical velocity field (for destructor)
  int    nSclr; //!< Number of scalar quantities to integrate
  int    maxit; //!< Maximum number of AL-iterations
  double eps[3]; //!< Convergence tolerances
};

#endif
