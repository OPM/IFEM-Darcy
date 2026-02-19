// $Id$
//==============================================================================
//!
//! \file main_DarcyCorrection.C
//!
//! \date Jan 20 2026
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for Darcy transport correction.
//!
//==============================================================================

#include "ASMmxBase.h"
#include "DarcyTransportCorr.h"
#include "DarcyArgs.h"
#include "SIMDarcyTransportCorr.h"

#include "ASMenums.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Profiler.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMoptions.h"
#include "SIMSolver.h"

#include <iostream>
#include <string>


/*!
  \brief Launch a simulator using a specified solver template.
  \param infile The input file to parse
  \param args Darcy arguments
*/

template<class Dim>
int runSimulator(char* infile, const DarcyArgs& args)
{
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS2;
  DarcyTransportCorr itg(Dim::dimension);
  const auto fields = args.mixed ? SIMinput::CharVec{Dim::dimension,2}
                                 : SIMinput::CharVec{Dim::dimension};
  SIMDarcyTransportCorr<Dim> darcy(itg, fields);
  SIMSolverStat solver(darcy);

  utl::profiler->start("Model input");

  if (!darcy.read(infile) || !solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (!darcy.preprocess())
    return 2;

  darcy.init();

  if (darcy.opt.dumpHDF5(infile))
    solver.handleDataOutput(darcy.opt.hdf5,darcy.getProcessAdm());

  int res = solver.solveProblem(infile,"Solving Darcy transport correction problem");
  // if (!res)
  //   darcy.printFinalNorms(TimeStep());

  return res;
}


/*!
  \brief Main program for the isogeometric Darcy solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  char* infile = nullptr;
  DarcyArgs args;

  IFEM::Init(argc,argv,"Darcy transport correction solver");
  for (int i = 1; i < argc; i++)
    if (args.parseArgComplex(argc,argv,i))
      ;
    else if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile) {
      infile = argv[i];
      if (!args.readXML(infile,false))
        return 1;
      i = 0;
    }
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>] [-hdf5]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu>] [-nv <nv>]"
              <<" [-nw <nw>]]\n";
    return 0;
  }

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;

  utl::profiler->stop("Initialization");

  if (args.dim == 3)
    return runSimulator<SIM3D>(infile,args);
  else if (args.dim == 2)
    return runSimulator<SIM2D>(infile,args);
  else
    return runSimulator<SIM1D>(infile,args);
}
