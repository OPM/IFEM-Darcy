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
#include "DarcyArgs.h"
#include "Profiler.h"


template<class Dim, template<class T> class Solver>
int runSimulator(char* infile)
{
  SIMDarcy<Dim> darcy;
  Solver<SIMDarcy<Dim>> solver(darcy);

  utl::profiler->start("Model input");

  if (!darcy.read(infile) || !solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (!darcy.preprocess())
    return 2;

  if (darcy.opt.dumpHDF5(infile))
    solver.handleDataOutput(darcy.opt.hdf5);

  return solver.solveProblem(infile,"Solving Darcy problem");
}


template<class Dim>
int runSimulator1(char* infile, bool adaptive)
{
  if (adaptive)
    return runSimulator<Dim, SIMSolverAdap>(infile);

  return runSimulator<Dim, SIMSolverStat>(infile);
}


int main(int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  char* infile = nullptr;
  DarcyArgs args;

  IFEM::Init(argc,argv,"Darcy solver");

  int ignoreArg = -1;
  for (int i = 1; i < argc; i++)
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
              <<"       [-ignore <p1> <p2> ...] [-fixDup]\n";
    return 0;
  }

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;
  utl::profiler->stop("Initialization");

  if (args.dim == 3)
    return runSimulator1<SIM3D>(infile,args.adap);
  else if (args.dim == 2)
    return runSimulator1<SIM2D>(infile,args.adap);
  else
    return runSimulator1<SIM1D>(infile,args.adap);
}
