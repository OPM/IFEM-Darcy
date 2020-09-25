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
#include "SIMargsBase.h"
#include "Profiler.h"


/*!
  \brief Launch a simulator using a specified solver template.
  \param infile The input file to parse
*/

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
    solver.handleDataOutput(darcy.opt.hdf5,darcy.getProcessAdm());

  return solver.solveProblem(infile,"Solving Darcy problem");
}

/*!
  \brief Choose a solver template and then launch a simulator.
  \param infile The input file to parse
  \param adaptive \e true to do an adaptive simulation
*/

template<class Dim>
int runSimulator1(char* infile, bool adaptive)
{
  if (adaptive)
    return runSimulator<Dim, SIMSolverAdap>(infile);

  return runSimulator<Dim, SIMSolverStat>(infile);
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
  SIMargsBase args("darcy");

  IFEM::Init(argc,argv,"Darcy solver");
  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
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
    return runSimulator1<SIM3D>(infile,args.adap);
  else if (args.dim == 2)
    return runSimulator1<SIM2D>(infile,args.adap);
  else
    return runSimulator1<SIM1D>(infile,args.adap);
}
