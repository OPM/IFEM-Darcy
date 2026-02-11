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
#include "SIMconfigure.h"

#include <string>

class DarcyTransportCorr;
class DataExporter;
class TimeStep;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Driver class for Darcy transport correction problems.
*/

template<class Dim>
class SIMDarcyTransportCorr : public Dim
{
public:
  //! \brief Setup properties.
  struct SetupProps {
    DarcyTransportCorr* itg = nullptr; //!< Pointer to integrand
  };

  //! \brief Default constructor.
  //! \param itg Integrand to use
  explicit SIMDarcyTransportCorr(DarcyTransportCorr& itg);

  //! \brief Mixed constructor.
  //! \param itg Integrand to use
  //! \param nf Fields on each basis
  explicit SIMDarcyTransportCorr(DarcyTransportCorr& itg,
                                 const std::vector<unsigned char>& nf);

  //! \brief Construct from setup properties.
  //! \param props Setup properties
  explicit SIMDarcyTransportCorr(const SetupProps& props)
    : SIMDarcyTransportCorr(*props.itg)
  {}

  //! \brief Destructor.
  ~SIMDarcyTransportCorr() override;

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "DarcyTransportCorr"; }

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

  //! \brief Computes the solution for the current time step.
  bool solveStep(const TimeStep& tp);

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solvec The solution vector
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  void printSolutionSummary(const Vector& solvec, int, const char*,
                            std::streamsize outPrec) override;

  //! \brief Adds a global multiplier for the constraint.
  //! \param nnod Number of nodes in model
  bool preprocessBeforeAsmInit(int& nnod) override;

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the generic Neumann property codes.
  void preprocessA() override;

protected:
  bool constrainIntegratedLag = false; //!< Constrain integrated multiplier
  Vector qSol; //!< Solution vector
  int vCode = 0; //!< Velocity anasol code
};


//! \brief Partial specialization for configurator.
template<class Dim>
struct SolverConfigurator<SIMDarcyTransportCorr<Dim>> {
  //! \brief Configure a SIMDarcyAdvection instance.
  //! \param ad The SIMDarcy instance to configure
  //! \param[in] props Configuration properties
  //! \param[in] infile The input file to read
  int setup(SIMDarcyTransportCorr<Dim>& ad,
            const typename SIMDarcyTransportCorr<Dim>::SetupProps& props,
            char* infile);
};

#endif
