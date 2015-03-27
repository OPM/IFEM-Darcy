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
#include "AdaptiveSIM.h"
#include "Utilities.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"

template<class Dim>
int runSimulator(const char* infile, bool adaptive)
{
  SIMDarcy<Dim> darcy;

  SIMinput* model = &darcy;
  AdaptiveSIM* aSim=NULL;
  if (adaptive)
    model = aSim = new AdaptiveSIM(&darcy);

  // read input file
  if (!model->read(infile))
    return 1;

  // configure finite element library
  if (!darcy.preprocess())
    return 2;

  if (aSim)
    aSim->setupProjections();

  // setup integration
  darcy.setQuadratureRule(model->opt.nGauss[0], true);
  darcy.initSystem(model->opt.solver, 1, 1);
  darcy.setAssociatedRHS(0, 0);
  darcy.setMode(SIM::STATIC);
  if (adaptive)
    aSim->initAdaptor(0,2);

  Vector sol;

  // HDF5 output
  DataExporter* exporter=NULL;
  if (model->opt.dumpHDF5(infile))
  {
    exporter = new DataExporter(true);
    int results = DataExporter::PRIMARY | DataExporter::SECONDARY;
    exporter->registerField("p,v", "primary", DataExporter::SIM, results);
    exporter->setFieldValue("p,v", &darcy, adaptive?&aSim->getSolution():&sol);
    exporter->registerWriter(new HDF5Writer(darcy.opt.hdf5,darcy.getProcessAdm()));
    exporter->registerWriter(new XMLWriter(darcy.opt.hdf5,darcy.getProcessAdm()));
    if (adaptive)
      exporter->setNormPrefixes(aSim->getNormPrefixes());
  }

  // VTF output

  bool iterate=true;
  for(int iStep=1;iterate;++iStep) {
    if (adaptive) {
      if (!aSim->solveStep(infile,iStep))
        return 5;
      else if (!aSim->writeGlv(infile,iStep,2))
        return 6;
    } else {
      // assemble the linear system
      if (!darcy.assembleSystem())
        return 3;

      // solve the linear system
      if (!darcy.solveSystem(sol, 1))
        return 4;

      // Evaluate solution norms
      Matrix eNorm;
      Vectors gNorm;
      darcy.setQuadratureRule(model->opt.nGauss[1]);
      if (!darcy.solutionNorms(Vectors(1,sol),Vectors(),eNorm,gNorm))
        return 4;

      darcy.printNorms(gNorm,IFEM::cout);

      // VTF output
      if (darcy.opt.format >= 0)
      {
        int geoBlk = 0, nBlock=0;
        darcy.writeGlvG(geoBlk, infile);
        darcy.writeGlvS(sol, 1, nBlock, 0.0, 0);
        //darcy.savePoints("")
        darcy.writeGlvStep(1, 0.0, 1);
        darcy.closeGlv();
      }
    }

    if (exporter)
      exporter->dumpTimeLevel(NULL, adaptive);
    if (adaptive)
      iterate = aSim->adaptMesh(iStep+1);
    else
      iterate = false;
  }

  delete exporter;
  return 0;
}


int main(int argc, char ** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  int  i;
  char ndim = 3;
  char* infile = 0;
  bool adaptive = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcmp(argv[i],"-adap"))
      adaptive = true;
    else if (!infile)
      infile = argv[i];
    else
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

  IFEM::cout <<"\n >>> IFEM Darcy equation solver <<<"
             <<"\n ====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator<SIM3D>(infile,adaptive);
  else if (ndim == 2)
    return runSimulator<SIM2D>(infile,adaptive);
  else
    return runSimulator<SIM1D>(infile,adaptive);

  return 1;
}
