/**
 * @file simulator.cpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Feb 2017
 * @section LICENSE
 * 
 * The GNU LGPL v3 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * pSSAlib Command Line Interface: simulator component
 */

#include <iostream>

#include <sbml/SBMLTypes.h>

#include "PSSA.h"
#include "util/FileSystem.h"

#include "CmdLineOptions.hpp"
using namespace pssalib::program_options;

// Debug backtrace
#ifdef __linux__
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

void handler(int sig)
{
  static const size_t szBacktrace = 10;
  
  void *array[szBacktrace];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, szBacktrace);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO); // write trace to stderr
  exit(1);
}
#endif

// Auxiliary struct
struct convertResults : public std::unary_function<UINTEGER, UINTEGER>
{
  UINTEGER operator()(const UINTEGER sr) const
  {
    UINTEGER result = 0;
    if(sr & ProgramOptionsSimulator::srTrajectory)
      result |= pssalib::datamodel::SimulationInfo::ofTrajectory;
    if(sr & ProgramOptionsSimulator::srFinalPopulations)
      result |= pssalib::datamodel::SimulationInfo::ofFinalPops;
    if(sr & ProgramOptionsSimulator::srTimePoints)
      result |= pssalib::datamodel::SimulationInfo::ofTimePoints;
    if(sr & ProgramOptionsSimulator::srTiming)
      result |= pssalib::datamodel::SimulationInfo::ofTiming;
    return result;
  }
};

/**
 * Main
 */
int main(int argc, char** argv)
{
#ifdef __linux__
  signal(SIGSEGV, handler);
#endif
  PSSALIB_MPI_IO_INIT;

  PSSALIB_MPI_COUT_OR_NULL << "\npSSAlib command line interface\n\n\n";

  // construct the SimulationInfo object
  pssalib::datamodel::SimulationInfo simInfo;

  // Parse command line and configuration file, if specified
  prog_opt::variables_map vm;
  ProgramOptionsSimulator poSimulator;
  if(!poSimulator.processCmdLineArgs(argc, argv, vm)) {
    return -127;
  }

  // Set up the simulation info and the pSSAlib instance.
  boost::scoped_ptr<pssalib::PSSA> ptrPSSA(new pssalib::PSSA());
  if (NULL == ptrPSSA.get()) {
    PSSALIB_MPI_CERR_OR_NULL << "Error : failed to allocate the simulation engine.\n";
    return -126;
  }

  simInfo.unSamplesTotal = poSimulator.getNumSamples();
  simInfo.dTimeStart = poSimulator.getTimeBegin();
  simInfo.dTimeStep = poSimulator.getTimeStep();
  simInfo.dTimeEnd = poSimulator.getTimeEnd();

  {
    convertResults cr;
    simInfo.unOutputFlags |= cr(poSimulator.getResultOutputFlags());
  }

  simInfo.eInitialPopulation = poSimulator.getInitialPopulation();
  simInfo.eBoundaryConditions = poSimulator.getBoundaryConditions();

  if(poSimulator.isSpatialSet())
  {
    std::pair<UINTEGER, UINTEGER> spatial = poSimulator.getSpatial();
    UINTEGER space[3] = { 0 };
    std::fill_n(space, spatial.first, spatial.second);
    simInfo.setDims(spatial.first, space);
  }

  if(poSimulator.isVerboseSet()&&poSimulator.isQuietSet()) {
    PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and "
      "'quiet' flags set, however, the later has priority over the former.\n";
  }

  // default output flags
  simInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog;

  if(poSimulator.isVerboseSet())
    simInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace;

  if(poSimulator.isQuietSet())
  {
    simInfo.unOutputFlags &= ~(pssalib::datamodel::SimulationInfo::ofTrace |
      pssalib::datamodel::SimulationInfo::ofInfo |
      pssalib::datamodel::SimulationInfo::ofWarning);

    if(poSimulator.getQuietCount() > 1)
    {
      simInfo.unOutputFlags &= ~(pssalib::datamodel::SimulationInfo::ofLog);
    }
    if(poSimulator.getQuietCount() > 2)
    {
      simInfo.unOutputFlags &= ~(pssalib::datamodel::SimulationInfo::ofStatus);
    }
    if(poSimulator.getQuietCount() > 3)
    {
      simInfo.unOutputFlags &= ~(pssalib::datamodel::SimulationInfo::ofError);
    }
  }
  else if(!poSimulator.isLogSet())
  {
    simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
  }

  if(poSimulator.isBenchmarkSet())
  {
      simInfo.unOutputFlags = ~(
        pssalib::datamodel::SimulationInfo::ofTrace      |
        pssalib::datamodel::SimulationInfo::ofInfo       |
        pssalib::datamodel::SimulationInfo::ofWarning    |
        pssalib::datamodel::SimulationInfo::ofError      |
        pssalib::datamodel::SimulationInfo::ofTrajectory |
        pssalib::datamodel::SimulationInfo::ofFinalPops  |
        pssalib::datamodel::SimulationInfo::ofTimePoints);
      simInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTiming;
      PSSALIB_MPI_COUT_OR_NULL << "Benchmarking, disable most outputs.\n";
  }

  // Simulator output
  if(!poSimulator.isQuietSet())
  {
    PSSALIB_MPI_COUT_OR_NULL << "Simulate model defined in :\n\tinput file : " 
                     << poSimulator.getInputFile() << "\nand output results to\n"
                     << "\tsave path : " << poSimulator.getOutputPath()
                     << "\nusing each of these methods:\n";
    std::transform(poSimulator.getMethods().begin(), poSimulator.getMethods().end(),
                   std::ostream_iterator<STRING>(PSSALIB_MPI_COUT_OR_NULL, "\t"), &pssalib::PSSA::getMethodName);
    PSSALIB_MPI_COUT_OR_NULL << std::endl;

    PSSALIB_MPI_COUT_OR_NULL << "Run " << simInfo.unSamplesTotal << " trials for '" << simInfo.dTimeEnd << "' seconds "
                     "outputting population every '" << simInfo.dTimeStep << "' seconds, beginning at '" 
                     << simInfo.dTimeStart << "' second(s).\n";

    PSSALIB_MPI_COUT_OR_NULL << "Output following resutls: \n";
    poSimulator.outputResultFlags(PSSALIB_MPI_COUT_OR_NULL,"\t");
    PSSALIB_MPI_COUT_OR_NULL << std::endl;

    if(NULL != poSimulator.getSpecies())
    {
      simInfo.pArSpeciesIds = new std::vector<STRING>(poSimulator.getSpecies()->size());
      std::copy(poSimulator.getSpecies()->begin(), poSimulator.getSpecies()->end(),
                simInfo.pArSpeciesIds->begin());
      if(simInfo.pArSpeciesIds->size() > 0)
      {
        PSSALIB_MPI_COUT_OR_NULL << "Output results only for species with following ids :\n";
        std::copy(simInfo.pArSpeciesIds->begin(), simInfo.pArSpeciesIds->end(),
                  std::ostream_iterator<STRING>(PSSALIB_MPI_COUT_OR_NULL, "\t"));
        PSSALIB_MPI_COUT_OR_NULL << std::endl;
      }
      else
      {
        PSSALIB_MPI_COUT_OR_NULL << "Do not output results for any species.\n";
      }
    }
    else
    {
      PSSALIB_MPI_COUT_OR_NULL << "Output results for all species.\n";
    }

    PSSALIB_MPI_COUT_OR_NULL << std::endl;
  }

  // parse the SBML model
  if (!simInfo.readSBMLFile(poSimulator.getInputFile()))
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error : failed to load SBML model from file '" << poSimulator.getInputFile() << "'.\n";
    return -125;
  }

  // create output directory
  STRING outputPathCurr;
  std::vector< STRING > arOutputPath;
  arOutputPath.push_back(poSimulator.getOutputPath());
  if(!pssalib::util::makeDir(arOutputPath,outputPathCurr))
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error : failed to create output directory '" << outputPathCurr << "'.\n";
    return -124;
  }

  {
    STRING cfgPath;
    pssalib::util::makeFilePath(poSimulator.getOutputPath(), STRING("pssa.cfg"), cfgPath);
    if(
#ifdef HAVE_MPI
      pssalib::getMPIWrapperInstance().isMaster() &&
#endif
      !serializeVariableMap(vm, cfgPath))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to store configuration to file '" << cfgPath << "'.\n";
      return -123;
    }
  }

  // Redirect output for species IDs
  FILESTREAMBUFFER fsbSpecies;
#ifdef HAVE_MPI
  if(pssalib::getMPIWrapperInstance().isMaster())
#endif
  {
    STRING idsPath;
    pssalib::util::makeFilePath(poSimulator.getOutputPath(),
                                pssalib::datamodel::SimulationInfo::arFileNames[
                                pssalib::datamodel::SimulationInfo::outputFlagToStreamIndex(pssalib::datamodel::SimulationInfo::ofSpeciesIDs)],
                                idsPath);

    if(NULL == fsbSpecies.open(idsPath.c_str(), std::ios_base::out | std::ios_base::trunc))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to open species ids file '" << idsPath << "'.\n";
      return -122;
    }
    simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofSpeciesIDs, &fsbSpecies);
  }

  ///////////////////////////////////////////
  // Run the simulations looping through SSAs
  bool bResult = true;
  for(std::vector<pssalib::PSSA::EMethod>::const_iterator mi = poSimulator.getMethods().begin();
      mi != poSimulator.getMethods().end(); mi++)
  {
    if(!ptrPSSA->setMethod((*mi)))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to set simulation method "
        << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi) << std::endl;
      continue;
    }

    // Create the output directory
    arOutputPath.clear();
    arOutputPath.push_back(poSimulator.getOutputPath());
    arOutputPath.push_back(pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi));
    if(!pssalib::util::makeDir(arOutputPath,outputPathCurr))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to create output directory '" << outputPathCurr << "'.\n";
      return -121;
    }
    simInfo.strOutput = outputPathCurr;

    if(!poSimulator.isQuietSet())
    {
      PSSALIB_MPI_COUT_OR_NULL << "simulating '" << simInfo.getModel().getName() << "' using "
                        << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi)
                        << "  ... \n";
    }

    // Run the simulation
    bool result = ptrPSSA->run(&simInfo);
    if(result)
    {
      if(!poSimulator.isQuietSet())
      {
        PSSALIB_MPI_COUT_OR_NULL << "done!\n";
      }
    }
    else
    {
      if(poSimulator.isQuietSet())
        std::cerr << "FAILED!\n";
      else
        std::cerr << "FAILED to simulate '" << simInfo.getModel().getName() << "' using "
                  << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi)
                  << "  ... \n";
    }

    bResult &= result;

    simInfo.unOutputFlags &= ~(pssalib::datamodel::SimulationInfo::ofSpeciesIDs);
  }

  if(bResult)
    return 0;
  else
    return -1;
}
