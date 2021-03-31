/**
 * @file CmdLineOptions.h
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
 * Auxiliary classes for parsing command line options
 */

#include <limits>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "util/ProgramOptionsBase.hpp"

#include <boost/format.hpp>

#include "PSSA.h"
#include "datamodel/detail/VolumeDecomposition.hpp"

#include "util/ProgramOptionsBase.hpp"

#ifndef PSSALIB_CLI_CMD_LINE_OPTIONS_H_
#define PSSALIB_CLI_CMD_LINE_OPTIONS_H_

namespace prog_opt = boost::program_options;
namespace prop_tree = boost::property_tree;

namespace pssalib
{
namespace program_options
{

//////////////////////////////
// Type for counting number of occurences of a switch on the command line
// Adopted with modification from http://stackoverflow.com/questions/31696328/boost-program-options-using-zero-parameter-options-multiple-times
typedef struct tagCLIOptionCounter
{
  UINTEGER count;

  tagCLIOptionCounter()
    : count(0)
  {
    // Do nothing
  };
} CLIOptionCounter;

void validate(boost::any& v, std::vector<std::string> const& xs, CLIOptionCounter *, long)
{
  try
  {
    if (v.empty()) v = CLIOptionCounter();
    ++boost::any_cast<CLIOptionCounter &>(v).count;
  }
  catch(const boost::bad_any_cast & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "CLIOptionCounter validator:" << e.what() << std::endl;
  }
}
//////////////////////////////

/**
 * Program options wrapper covering options common for both
 * simulator and analyse components of the pSSAlib CLI.
 */
template< bool RequireMethods >
class ProgramOptionsCommon : public ProgramOptionsBase
{
//////////////////////////////
// Attributes
protected:
  //! Quiet or verbose?
  UINTEGER m_unQuietCount;
  bool     m_bVerbose;

  //! Output path
  STRING m_strOutputPath;

  //! Species ids to process
  std::vector<STRING> *m_pArSpeciesIds;

  //! Time point when the output begins
//   bool m_bIsTimeBeginSet;
  REAL m_dTimeBegin;

  //! Final time point of the simulation
//   bool m_bIsTimeEndSet;
  REAL m_dTimeEnd;

  //! Number of trajectories
  UINTEGER m_unSamples;

  //! Simulation methods
  std::vector<pssalib::PSSA::EMethod> m_arMethods;

//////////////////////////////
// Constructors
public:

  //! Constructor
  ProgramOptionsCommon(const char *s)
    : ProgramOptionsBase(s)
    , m_pArSpeciesIds(NULL)
  {
    resetToDefault();
  }

  //! Destructor
  ~ProgramOptionsCommon()
  {
    if(NULL != m_pArSpeciesIds)
      delete m_pArSpeciesIds;
  }

//////////////////////////////
// Methods
protected:

  /**
   * Initiliaze the description object.
   * @return @true if initialization is successfull, @false otherwise.
   */
  virtual bool initialize()
  {
    if(!ProgramOptionsBase::initialize())
      return false;

    try
    {
      m_poDesc.add_options()
        ("output-path,o",   prog_opt::value<STRING>()->required(),                  "Output path")
        ("species,s",       prog_opt::value< CLIOptionCommaSeparatedList >(),       "Comma-separated list of species ids for which output "
                                                                                    "is desired. Use an empty string if none. If omitted, "
                                                                                    "all species are output.")
        ("tend",            prog_opt::value<REAL>(),                                "End time of the simulation")
        ("tstart",          prog_opt::value<REAL>()->default_value(0.0),            "Time point when output of trajectories begins")
        ("num-samples,n",   prog_opt::value<UINTEGER>(),                            "Number of samples to collect")
        ("methods,m",       RequireMethods ?
                            prog_opt::value< CLIOptionCommaSeparatedList >()->required() :
                            prog_opt::value< CLIOptionCommaSeparatedList >(),       "A comma-separated list of simulation method ids:"
                                                                                    "\n0,dm - Gillespie's Direct Method"
                                                                                    "\n1,pdm - Partial Propensity Direct Method"
                                                                                    "\n2,pssacr - pSSA with Composition-Rejection Sampling"
                                                                                    "\n3,spdm - Sorting Partial Propensity Direct Method")
        ("verbose,v",                                                               "Output additional information about the simulation")
        ("quiet,q",         prog_opt::value<CLIOptionCounter>()->zero_tokens(),     "Output only essential information about the simulation, "
                                                                                    "may be specified multiple times for a cumulative effect")
        ;

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

//////////////////////////////
// Methods
public:

  /**
   * Reset the object with pre-set values.
   */
  virtual void resetToDefault()
  {
    m_unQuietCount = 0;
    m_bVerbose = false;

    m_strOutputPath.clear();

    if(NULL != m_pArSpeciesIds)
      delete m_pArSpeciesIds;
    m_pArSpeciesIds = NULL;

    m_dTimeBegin = -1.0 * std::numeric_limits<REAL>::infinity(); // < 0 => not set
    m_dTimeEnd = -1.0 * std::numeric_limits<REAL>::infinity(); // < 0 => not set

    m_unSamples = 0;

    m_arMethods.clear();
  }

  virtual bool parseVariableMap(prog_opt::variables_map & vm)
  {
    typedef std::map<STRING, UINTEGER> MAPPING_TYPE;
    typedef std::vector<UINTEGER> RESULT_TYPE;

    MAPPING_TYPE mapping;
    RESULT_TYPE result;

    try
    {
      if(0 == vm.count("output-path"))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : no output path specified." << std::endl;
        return false;
      }
      else
        m_strOutputPath = vm["output-path"].as<STRING>();

      if(vm.count("species") > 0)
      {
        if(NULL != m_pArSpeciesIds)
          delete m_pArSpeciesIds;
        m_pArSpeciesIds = new std::vector<STRING>();

        CLIOptionCommaSeparatedList species = vm["species"].as< CLIOptionCommaSeparatedList >();
        species.parse(*m_pArSpeciesIds, true, false, false);
      }
      else // default
        m_pArSpeciesIds = NULL;

      if(0 == vm.count("methods"))
      {
        if(RequireMethods)
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error : no methods specified." << std::endl;
          return false;
        }
      }
      else
      {
        mapping.clear();
        result.clear();

        mapping[STRING("0")] = pssalib::PSSA::M_DM;
        mapping[STRING("dm")] = pssalib::PSSA::M_DM;
        mapping[STRING("1")] = pssalib::PSSA::M_PDM;
        mapping[STRING("pdm")] = pssalib::PSSA::M_PDM;
        mapping[STRING("2")] = pssalib::PSSA::M_PSSACR;
        mapping[STRING("pssacr")] = pssalib::PSSA::M_PSSACR;
        mapping[STRING("3")] = pssalib::PSSA::M_SPDM;
        mapping[STRING("spdm")] = pssalib::PSSA::M_SPDM;

        CLIOptionCommaSeparatedList methods = vm["methods"].as< CLIOptionCommaSeparatedList >();
        methods.parse(mapping, result, true, false, false);

        if(0 == result.size())
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid method specification. Valid values are:\n\n";
          std::for_each(mapping.begin(), mapping.end(),
                        printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
          PSSALIB_MPI_CERR_OR_NULL << "\n\n";
          return false;
        }
        else
        {
          m_arMethods.resize(result.size());
          std::transform(result.begin(), result.end(), m_arMethods.begin(),
                         converTypes<UINTEGER, pssalib::PSSA::EMethod>());
        }
      }

      if(vm.count("quiet") > 0)
      {
        CLIOptionCounter quietCounter = vm["quiet"].as<CLIOptionCounter>();
        m_unQuietCount = quietCounter.count;
      }
      m_bVerbose = (vm.count("verbose") > 0);

      if(vm.count("tstart") > 0)
        m_dTimeBegin = std::max(vm["tstart"].as<REAL>(), 0.0);
      if(vm.count("tend") > 0)
        m_dTimeEnd = std::max(vm["tend"].as<REAL>(), 0.0); // (vm.count("tstart") > 0) ? m_dTimeBegin : 

      if((vm.count("tstart") > 0)&&(vm.count("tend") > 0))
      {
        if(getTimeBegin() > getTimeEnd())
          throw prog_opt::error("intial time point cannot lie past the final one!");

        if(getTimeBegin() == getTimeEnd())
          throw prog_opt::error("intial and final time points cannot coinside!");
      }

      if(vm.count("num-samples") > 0)
        m_unSamples = vm["num-samples"].as<UINTEGER>();

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
    catch (boost::bad_any_cast &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

  //////////////////////////////
  // Getters & setters
  //////////////////////////////

  bool isVerboseSet() const
  {
    return m_bVerbose;
  }

  bool isQuietSet() const
  {
    return (0 != m_unQuietCount);
  }

  UINTEGER getQuietCount() const
  {
    return m_unQuietCount;
  }

  const STRING & getOutputPath() const
  {
    return m_strOutputPath;
  }

  bool isTimeBeginSet() const
  {
    return !isinf(m_dTimeBegin);
  }

  REAL getTimeBegin() const
  {
    if(isinf(m_dTimeBegin))
      return 0.0;
    return m_dTimeBegin;
  }

  bool isTimeEndSet() const
  {
    return !isinf(m_dTimeEnd);
  }

  REAL getTimeEnd() const
  {
    if(isinf(m_dTimeEnd))
      return 0.0;
    return m_dTimeEnd;
  }

  bool isNumSamplesSet() const
  {
    return (0 != m_unSamples);
  }

  UINTEGER getNumSamples() const
  {
    if(0 == m_unSamples)
      return 10;
    return m_unSamples;
  }

  const std::vector<pssalib::PSSA::EMethod> & getMethods() const
  {
    return m_arMethods;
  }

  const std::vector<STRING> * getSpecies() const
  {
    return m_pArSpeciesIds;
  }
};

/**
 * Program options wrapper for the pSSAlib CLI simulator.
 */
class ProgramOptionsSimulator : public ProgramOptionsCommon<true>
{
//////////////////////////////
// Data type
public:

  typedef enum tagSimulatorResults
  {
    srTrajectory       = 0x01,
    srFinalPopulations = 0x02,
    srTimePoints       = 0x04,
    srTiming           = 0x08
  } SimulatorResults;

//////////////////////////////
// Attributes
protected:
  //! Log the simulation engine output to a file
  bool m_bLog;

  //! Is this a benchmark run?
  bool m_bBenchmark;

  //! Input SBML file
  STRING m_strInputFile;

  //! Time step
  bool m_bIsTimeStepSet;
  REAL m_dTimeStep;

  //! Volume override
  bool m_bIsTotalVolumeSet;
  REAL m_dTotalVolume;

  //! Initial population type
  //bool m_bIsTotalVolumeSet;
  pssalib::datamodel::detail::InitialPopulationType
    m_InitPop;

  //! Boundary conditions
  //bool m_bIsTotalVolumeSet;
  pssalib::datamodel::detail::BoundaryConditionsType
    m_BndCond;

  //! Spatial dimensions
  UINTEGER m_unSubreactors,
           m_unDimensions;

  //! Results to output
  UINTEGER m_unResults;

//////////////////////////////
// Constructors
public:

  //! Constructor
  ProgramOptionsSimulator()
    : ProgramOptionsCommon<true>("Generic simulator options")
  {
    resetToDefault();
  }

  //! Destructor
  ~ProgramOptionsSimulator()
  {
    // Do nothing
  }

//////////////////////////////
// Methods
protected:

  /**
   * Initiliaze the description object.
   * @return @true if initialization is successfull, @false otherwise.
   */
  virtual bool initialize()
  {
    if(!ProgramOptionsCommon<true>::initialize())
      return false;

    try
    {
      m_poDesc.add_options()
        ("dt",              prog_opt::value<REAL>()->default_value(0.1),            "Time interval between outputs")
        ("sbml-file,i",     prog_opt::value<STRING>()->required(),                  "SBML input file")
        ("results,r",       prog_opt::value< CLIOptionCommaSeparatedList >(),       "A comma-separated list of results to output:"
                                                                                    "\n0,trajectory - Trajectory of species population"
                                                                                    "\n1,finalVals - Populations at final time (used to compute pdfs)"
                                                                                    "\n2,timePoints - Output the time points to a separate file"
                                                                                    "\n3,timing - Output timing info (only useful if benchmarking is on)")
        ("total-volume",    prog_opt::value<REAL>()->default_value(1.0),            "Size of the total volume")
        ("bndcond",         prog_opt::value< CLIOptionCommaSeparatedList >(),       "Boundary conditions, can be either:"
                                                                                    "\n0,\"periodic\""
                                                                                    "\n1,\"reflexive\"")
        ("dimensions",      prog_opt::value< CLIOptionDimensionList >(),            "Number of subreactors followed by number of spatial dimensions separated "
                                                                                    "by an \"x\". The number  before \"x\" represents the number of subreactors "
                                                                                    "and the number following it is the number of spatial dimensions, can be , "
                                                                                    "e.g., 5x2 denotes a square lattice with 5 subreactors along each dimension.")
        ("initpop",         prog_opt::value< CLIOptionCommaSeparatedList >(),       "Specifies how the initial population specified in the SBML file is distributed across the volume:"
                                                                                    "\n0,\"distribute\" - the population is evenly distributed, i.e. each subvolume gets pop/num_volumes"
                                                                                    "\n1,\"concentrate\" - the population is concentrated in the middle cell, i.e. one subvolume gets everything"
                                                                                    "\n2,\"multiply\" - the population is multiplied, i.e. each subvolume gets the total population")
        ("log,l",                                                                   "Log simulation engine output to a file in the output subdir")
        ("benchmark,b",                                                             "Benchmark the algorithm (suppresses most outputs and produces timing data)")
        ;

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

//////////////////////////////
// Methods
public:

  /**
   * Convert result flags into strings & output them to a stream
   * @param os Output stream
   * @param delim Value delimiter
   */
  void outputResultFlags(std::ostream & os, const char * delim ="\t")
  {
    if(m_unResults & srTrajectory)
      os << "'Population Trajectories'" << delim;
    if(m_unResults & srFinalPopulations)
      os << "'Populations at Final Timepoints'" << delim;
    if(m_unResults & srTimePoints)
      os << "'Time Points'" << delim;
    if(m_unResults & srTiming)
      os << "'Timing'" << delim;
  }

  /**
   * Reset the object with pre-set values.
   */
  virtual void resetToDefault()
  {
    ProgramOptionsCommon<true>::resetToDefault();

    m_dTimeStep = std::numeric_limits<REAL>::min(); // < 0 => not set

    m_bLog = m_bBenchmark = false;

    m_dTotalVolume = std::numeric_limits<REAL>::min(); // < 0 => not set

    m_InitPop = pssalib::datamodel::detail::IP_Invalid;
    m_BndCond = pssalib::datamodel::detail::BC_Invalid;

    m_unResults = 0;
    m_unSubreactors = 0;
    m_unDimensions = 0;
  }

  virtual bool parseVariableMap(prog_opt::variables_map & vm)
  {
    if(!ProgramOptionsCommon<true>::parseVariableMap(vm))
      return false;

    typedef std::map<STRING, UINTEGER> MAPPING_TYPE;
    typedef std::vector<UINTEGER> RESULT_TYPE;

    MAPPING_TYPE mapping;
    RESULT_TYPE result;

    try
    {
      if(0 == vm.count("sbml-file"))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : no SBML model file specified." << std::endl;
        return false;
      }
      else
        m_strInputFile = vm["sbml-file"].as<STRING>();

      if(vm.count("results") > 0)
      {
        mapping.clear();
        result.clear();

        mapping[STRING("0")] = srTrajectory;//pssalib::datamodel::SimulationInfo::ofTrajectory;
        mapping[STRING("trajectories")] = srTrajectory;//pssalib::datamodel::SimulationInfo::ofTrajectory;
        mapping[STRING("1")] = srFinalPopulations;//pssalib::datamodel::SimulationInfo::ofFinalPops;
        mapping[STRING("finalVals")] = srFinalPopulations;//pssalib::datamodel::SimulationInfo::ofFinalPops;
        mapping[STRING("2")] = srTimePoints;//pssalib::datamodel::SimulationInfo::ofTimePoints;
        mapping[STRING("timePoints")] = srTimePoints;//pssalib::datamodel::SimulationInfo::ofTimePoints;
        mapping[STRING("3")] = srTiming;//pssalib::datamodel::SimulationInfo::ofTiming;
        mapping[STRING("timing")] = srTiming;//pssalib::datamodel::SimulationInfo::ofTiming;

        CLIOptionCommaSeparatedList results = vm["results"].as< CLIOptionCommaSeparatedList >();
        results.parse(mapping, result, true, true, false);

        if(0 == result.size())
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid results specifications. Valid values are:\n\n";
          std::for_each(mapping.begin(), mapping.end(),
                        printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
          PSSALIB_MPI_CERR_OR_NULL << "\n\n";
          return false;
        }
        else
        {
          for(RESULT_TYPE::iterator it = result.begin(); it != result.end(); it++)
            m_unResults |= (*(it));
        }
      }
      else // default
        m_unResults = srTrajectory;

      if(vm.count("bndcond") > 0)
      {
        mapping.clear();
        result.clear();

        mapping[STRING("0")] = pssalib::datamodel::detail::BC_Periodic;
        mapping[STRING("periodic")] = pssalib::datamodel::detail::BC_Periodic;
        mapping[STRING("1")] = pssalib::datamodel::detail::BC_Reflexive;
        mapping[STRING("reflexive")] = pssalib::datamodel::detail::BC_Reflexive;

        CLIOptionCommaSeparatedList boundaryConditions = vm["bndcond"].as< CLIOptionCommaSeparatedList >();
        boundaryConditions.parse(mapping, result, false, true, true);

        if(0 == result.size())
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid boundary conditions. Valid values are:\n\n";
          std::for_each(mapping.begin(), mapping.end(),
                        printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
          PSSALIB_MPI_CERR_OR_NULL << "\n\n";
          return false;
        }
        else
        {
          m_BndCond = (pssalib::datamodel::detail::BoundaryConditionsType)(*(result.begin()));
        }
      }
      else // default
        m_BndCond = pssalib::datamodel::detail::BC_Periodic;

      // Spatial dimensions
      if(vm.count("dimensions") > 0)
      {
        result.clear();
        
        CLIOptionDimensionList dims = vm["dimensions"].as< CLIOptionDimensionList >();
        dims.parse(result, false, true, true);
        
        if(result.size() < 2)
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid definition for spatial structure. "
                                      "Expected number of subreactors followed by number"
                                      " of spatial dimensions separated by an \"x\".\n";
          return false;
        }
        else
        {
          m_unSubreactors = result[0];
          m_unDimensions  = result[1];

          if(0 == m_unSubreactors)
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error: invalid definition for spatial structure. "
                                        "Number of subreactors must be positive.\n";
            return false;
          }

          if(0 == m_unDimensions)
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error: invalid definition for spatial structure. "
                                        "Number of spatial dimensions must be positive.\n";
            return false;
          }

          if(result.size() > 2)
          {
            PSSALIB_MPI_CERR_OR_NULL << "Warning: additional members in definition "
                                        "for spatial structure were ignored.\n";
            return false;
          }
        }
      }
      else // default
        m_unSubreactors = m_unDimensions = 0;

      if(vm.count("initpop") > 0)
      {
        mapping.clear();
        result.clear();

        mapping[STRING("0")] = pssalib::datamodel::detail::IP_Distribute;
        mapping[STRING("distribute")] = pssalib::datamodel::detail::IP_Distribute;
        mapping[STRING("1")] = pssalib::datamodel::detail::IP_Concentrate;
        mapping[STRING("concentrate")] = pssalib::datamodel::detail::IP_Concentrate;
        mapping[STRING("2")] = pssalib::datamodel::detail::IP_Multiply;
        mapping[STRING("multiply")] = pssalib::datamodel::detail::IP_Multiply;

        CLIOptionCommaSeparatedList initpop = vm["initpop"].as< CLIOptionCommaSeparatedList >();
        initpop.parse(mapping, result, false, true, true);

        if(0 == result.size())
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid intial population specification. Valid values are:\n\n";
          std::for_each(mapping.begin(), mapping.end(),
                        printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
          PSSALIB_MPI_CERR_OR_NULL << "\n\n";
          return false;
        }
        else
        {
          m_InitPop = (pssalib::datamodel::detail::InitialPopulationType)(*(result.begin()));
        }
      }
      else // default
        m_InitPop = pssalib::datamodel::detail::IP_Distribute;

      m_bLog = (vm.count("log") > 0);
      m_bBenchmark = (vm.count("benchmark") > 0);

      m_dTimeStep = vm["dt"].as<REAL>();

      if(vm.count("total-volume"))
        m_dTotalVolume = vm["total-volume"].as<REAL>();

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
    catch (boost::bad_any_cast &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

  //////////////////////////////
  // Getters & setters
  //////////////////////////////
  
  bool isLogSet() const
  {
    return m_bLog;
  }
  
  bool isBenchmarkSet() const
  {
    return m_bBenchmark;
  }
  
  const STRING & getInputFile() const
  {
    return m_strInputFile;
  }
  
  REAL getTimeStep() const
  {
    return m_dTimeStep;
  }
  
  REAL getTotalVolume() const
  {
    return m_dTotalVolume;
  }
  
  bool isSpatialSet() const
  {
    return (0 != m_unDimensions*m_unSubreactors);
  }
  
  std::pair<UINTEGER, UINTEGER> getSpatial() const
  {
    return std::pair<UINTEGER, UINTEGER>(m_unDimensions, m_unSubreactors);
  }
  
  pssalib::datamodel::detail::
  InitialPopulationType getInitialPopulation() const
  {
    return m_InitPop;
  }
  
  pssalib::datamodel::detail::
  BoundaryConditionsType getBoundaryConditions() const
  {
    return m_BndCond;
  }
  
  UINTEGER getResultOutputFlags() const
  {
    return m_unResults;
  }
};

/**
 * Program options wrapper for the pSSAlib CLI analyzer.
 */
class ProgramOptionsAnalyzer : public ProgramOptionsCommon<false>
{
//////////////////////////////
// Data type
public:

  typedef enum tagAnalyzerResults
  {
    arTrajectory        = 0x01,
    arAverageTrajectory = 0x02,
    arPDF               = 0x04,
    arTiming            = 0x08
  } AnalyzerResults;
  
  typedef enum tagAnalyzerFormats
  {
    afCSV     = 0x01,
    afGnuPlot = 0x02,
    afVTK     = 0x04
  } AnalyzerFormats;

//////////////////////////////
// Attributes
protected:
  //! Log the simulation engine output to a file
  bool m_bLog;

  //! Output path
  STRING m_strInputPath;

  //! Analyzer output flags
  UINTEGER m_unResults;

  //! Output format
  AnalyzerFormats m_Format;

//////////////////////////////
// Constructors
public:

  //! Constructor
  ProgramOptionsAnalyzer()
    : ProgramOptionsCommon<false>("Generic analyzer options")
  {
    resetToDefault();
  }

  //! Destructor
  ~ProgramOptionsAnalyzer()
  {
    // Do nothing
  }

//////////////////////////////
// Methods
protected:

  /**
   * Initiliaze the description object.
   * @return @true if initialization is successfull, @false otherwise.
   */
  virtual bool initialize()
  {
    if(!ProgramOptionsCommon<false>::initialize())
      return false;

    try
    {
      m_poDesc.add_options()
        ("input-path,i",   prog_opt::value<STRING>()->required(),                  "Output path")
        ("results,r",      prog_opt::value< CLIOptionCommaSeparatedList >(),       "A comma-separated list of results to output:"
                                                                                   "\n0,trajectory - Trajectory of species population"
                                                                                   "\n1,avg-trajectory - Average trajectory of species population"
                                                                                   "\n2,pdf - Compute the PDF using populations at final time"
                                                                                   "\n3,timing - Analyze timing information obtained from benchmarking")
        ("format,f",       prog_opt::value<STRING>()->default_value(STRING("csv")),"Output format, one of \"csv\", \"gnuplot\" or \"vtk\".")
//         ("log,l",                                                                  "Log simulation engine output to a file in the output subdir")
        ;

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

//////////////////////////////
// Methods
public:
  
  /**
   * Convert result flags into strings & output them to a stream
   * @param os Output stream
   * @param delim Value delimiter
   */
  void outputResultFlags(std::ostream & os, const char * delim ="\t")
  {
    if(m_unResults & arTrajectory)
      os << "'Population Trajectories'" << delim;
    if(m_unResults & arAverageTrajectory)
      os << "'Average Trajectory'" << delim;
    if(m_unResults & arPDF)
      os << "'Probability Distribution Function'" << delim;
    if(m_unResults & arTiming)
      os << "'Timing'" << delim;
  }

  /**
   * Reset the object with pre-set values.
   */
  virtual void resetToDefault()
  {
    ProgramOptionsCommon<false>::resetToDefault();

    m_bLog = false;

    m_unResults = 0;
  }

  virtual bool parseVariableMap(prog_opt::variables_map & vm)
  {
    if(!ProgramOptionsCommon<false>::parseVariableMap(vm))
      return false;

    typedef std::map<STRING, UINTEGER> MAPPING_TYPE;
    typedef std::vector<UINTEGER> RESULT_TYPE;

    MAPPING_TYPE mapping;
    RESULT_TYPE result;

    try
    {
      if(0 == vm.count("input-path"))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : no input path specified." << std::endl;
        return false;
      }
      else
        m_strInputPath = vm["input-path"].as<STRING>();

      if(vm.count("results") > 0)
      {
        mapping.clear();
        result.clear();

        mapping[STRING("0")] = arTrajectory;
        mapping[STRING("trajectories")] = arTrajectory;
        mapping[STRING("1")] = arAverageTrajectory;
        mapping[STRING("average-trajectory")] = arAverageTrajectory;
        mapping[STRING("2")] = arPDF;
        mapping[STRING("pdf")] = arPDF;;
        mapping[STRING("3")] = arTiming;
        mapping[STRING("timing")] = arTiming;

        CLIOptionCommaSeparatedList results = vm["results"].as< CLIOptionCommaSeparatedList >();
        results.parse(mapping, result, true, true, true);

        if(0 == result.size())
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid results specification. Valid values are:\n\n";
          std::for_each(mapping.begin(), mapping.end(),
                        printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
          PSSALIB_MPI_CERR_OR_NULL << "\n\n";
          return false;
        }
        else
        {
          for(RESULT_TYPE::iterator it = result.begin(); it != result.end(); it++)
            m_unResults |= (*(it));
        }
      }
      else // default
        m_unResults = pssalib::datamodel::SimulationInfo::ofTrajectory;

      if(vm.count("format") > 0)
      {
        STRING strFormat = vm["format"].as< STRING >();

        if(0 == strFormat.compare("csv"))
          m_Format = afCSV;
        else if(0 == strFormat.compare("gnuplot"))
           m_Format = afGnuPlot;
        else if(0 == strFormat.compare("vtk"))
           m_Format = afVTK;
        else
        {
          PSSALIB_MPI_CERR_OR_NULL << "Error: invalid format specification \"" << strFormat
                                   <<"\". Valid values are: \"csv\", \"gnuplot\" and "
                                   "\"vtk\".\n\n";
          return false;
        }
      }
      else // default
        m_Format = afCSV;

      return true;
    }
    catch (prog_opt::error &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
    catch (boost::bad_any_cast &e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
      return false;
    }
  }

  //////////////////////////////
  // Getters & setters
  //////////////////////////////

  bool isLogSet() const
  {
    return m_bLog;
  }

  AnalyzerFormats getFormat() const
  {
    return m_Format;
  }

  const STRING & getInputPath() const
  {
    return m_strInputPath;
  }

  UINTEGER getResultOutputFlags() const
  {
    return m_unResults;
  }
};

void printVariableMap(prog_opt::variables_map & vm, std::ostream & os)
{
  try
  {
    for (prog_opt::variables_map::iterator it = vm.begin(); it != vm.end(); it++ )
    {
      if(it->second.defaulted())
        continue;

      if(it->second.empty())
      {
        os << it->first << " => " << "true" << std::endl;
      }
      else if(typeid(UINTEGER) == it->second.value().type())
      {
        os << it->first << " => " << boost::any_cast<UINTEGER>(it->second.value()) << std::endl;
      }
      else if(typeid(INTEGER) == it->second.value().type())
      {
        os << it->first << " => " << boost::any_cast<INTEGER>(it->second.value()) << std::endl;
      }
      else if(typeid(REAL) == it->second.value().type())
      {
        os << it->first << " => " << boost::any_cast<REAL>(it->second.value()) << std::endl;
      }
      else if(typeid(CLIOptionCommaSeparatedList) == it->second.value().type())
      {
        os << it->first << " => "
           << (STRING)boost::any_cast< CLIOptionCommaSeparatedList >(it->second.value()) << std::endl;
      }
      else if(typeid(CLIOptionCounter) == it->second.value().type())
      {
        os << it->first << " => "
           << boost::any_cast<CLIOptionCounter>(it->second.value()).count << std::endl;
      }
      else
      {
        os << it->first << " => " << boost::any_cast<STRING>(it->second.value()) << std::endl;
      }
    }
  }
  catch(const boost::bad_any_cast &e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error printing variable map: " << e.what() << std::endl;
  }
}

bool serializeVariableMap(prog_opt::variables_map & vm, STRING & strFileName)
{
  bool bResult = true;

  {
    std::ofstream ofsIniFile;

    ofsIniFile.open(strFileName.c_str());
    if(!ofsIniFile.good())
    {
      ofsIniFile.close();
      PSSALIB_MPI_CERR_OR_NULL << "Error during serialisation: could not open file '"
        << strFileName << "'." << std::endl;
      bResult = false;
    }

    if(bResult)
    {
      try
      {
        prop_tree::ptree root;

        for (prog_opt::variables_map::iterator it = vm.begin(); it != vm.end(); it++ )
        {
          if(it->second.defaulted())
            continue;

          if(it->second.empty())
          {
            root.put(it->first, "true");
          }
          else if(typeid(UINTEGER) == it->second.value().type())
          {
            root.put(it->first, boost::format("%d") % boost::any_cast<UINTEGER>(it->second.value()));
          }
          else if(typeid(INTEGER) == it->second.value().type())
          {
            root.put(it->first, boost::format("%d") % boost::any_cast<INTEGER>(it->second.value()));
          }
          else if(typeid(REAL) == it->second.value().type())
          {
            root.put(it->first, boost::format("%f") % boost::any_cast<REAL>(it->second.value()));
          }
          else if(typeid(CLIOptionCommaSeparatedList) == it->second.value().type())
          {
            root.put(it->first, (STRING)boost::any_cast< CLIOptionCommaSeparatedList >(it->second.value()));
          }
          else if(typeid(CLIOptionCounter) == it->second.value().type())
          {
            root.put(it->first, boost::any_cast<CLIOptionCounter>(it->second.value()).count);
          }
          else
          {
            root.put(it->first, boost::any_cast<STRING>(it->second.value()));
          }
        }

        prop_tree::write_ini(ofsIniFile, root);
      }
      catch(const boost::bad_any_cast &e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error during serialisation: " << e.what() << std::endl;
        bResult = false;
      }
      catch(const prop_tree::ini_parser_error & e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error during serialisation: " << e.what() << std::endl;
        bResult = false;
      }
    }

    if(bResult&&!ofsIniFile.good())
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error during serialisation: could not open file '"
        << strFileName << "'." << std::endl;
      bResult = false;
    }

    ofsIniFile.close();
  }

  return bResult;
}

} } // close program_options & pssalib namespaces

#endif /* PSSALIB_CLI_CMD_LINE_OPTIONS_H_ */
