// $Id$
//==============================================================================
//!
//! \file main_Darcy.C
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Main program for an isogeometric solver for Darcy flow.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMDarcy.h"
#include "SIMSolverAdap.h"
#include "Utilities.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "AppCommon.h"
#include "TimeStep.h"
#include "DarcyArgs.h"


template<class Dim, template<class T> class Solver>
int runSimulator(char* infile)
{
  SIMDarcy<Dim> darcy;
  Solver<SIMDarcy<Dim>> solver(darcy);
  int res = ConfigureSIM(darcy, infile, false);

  if (res)
    return res;

  // HDF5 output
  std::unique_ptr<DataExporter> exporter;

  if (darcy.opt.dumpHDF5(infile))
    exporter.reset(SIM::handleDataOutput(darcy, solver, darcy.opt.hdf5));

  return solver.solveProblem(infile, exporter.get(), "Solving Darcy problem", false);
}


template<class Dim>
int runSimulator1(char* infile, bool adaptive)
{
  if (adaptive)
    return runSimulator<Dim, SIMSolverAdap>(infile);

  return runSimulator<Dim, SIMSolver>(infile);
}


int main(int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int  i;
  char* infile = 0;
  DarcyArgs args;

  IFEM::Init(argc,argv,"Darcy solver");

int ignoreArg = -1;
  for (i = 1; i < argc; i++)
    if (i == ignoreArg || SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      args.dim = 2;
    else if (!strcmp(argv[i],"-1D"))
      args.dim = 1;
    else if (!strcmp(argv[i],"-adap"))
      args.adap = true;
    else if (!infile) {
      infile = argv[i];
      ignoreArg = i;
      if (!args.readXML(infile,false))
        return 1;
      i = 0;
    } else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
              <<" [-free] [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>]"
              <<"\n       [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
              <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
              <<"       [-ignore <p1> <p2> ...] [-fixDup]" << std::endl;
    return 0;
  }

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (args.dim == 3)
    return runSimulator1<SIM3D>(infile,args.adap);
  else if (args.dim == 2)
    return runSimulator1<SIM2D>(infile,args.adap);
  else
    return runSimulator1<SIM1D>(infile,args.adap);
}
