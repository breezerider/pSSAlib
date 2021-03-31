/**
 * @file main.cpp
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
 * Benchmarks
 */

#include <iostream>

// #include <boost/tokenizer.hpp>
// #include <boost/lexical_cast.hpp>
// #include <boost/algorithm/string.hpp>

#include "PSSA.h"

#include "util/MPIWrapper.h"
#include "util/ProgramOptionsBase.hpp"

#include "CyclicLinearChain.hpp"
#include "ColloidalAggregation.hpp"

using namespace pssalib::program_options;

// enumerate all test case ids
enum tagTestCases
{
  tcCLC = 0x0001,
  tcCA  = 0x0002,
  tcAll = 0x0003
};

class Benchmarks : public ProgramOptionsBase
{
//////////////////////////////
// Attributes
private:

  //! Output options
  bool m_bVerbose, m_bQuiet;

  //! Simulation methods
  UINTEGER m_unMethods;

  //! Test cases
  UINTEGER m_unTests;

  //! End time
  REAL     m_dTimeEnd;

  //! Number of samples
  UINTEGER m_unSamples;

  // network sizes
  std::vector< UINTEGER > m_arunSizes;

//////////////////////////////
// Constructors
public:

  //! Constructor
  Benchmarks()
    : ProgramOptionsBase("Options for Benchmarks")
  {
    resetToDefault();
  }

  //! Destructor
  ~Benchmarks()
  {
    // Do nothing
  }

//////////////////////////////
// Methods
protected:

  /**
   * @copydoc ProgramOptionsBase::initialize()
   */
  virtual bool initialize()
  {
    // call base class member function
    if(!ProgramOptionsBase::initialize())
      return false;

    try
    {
      m_poDesc.add_options()
        ("tend",              prog_opt::value<REAL>()->default_value(1000.0),         "End time of the simulation")
        ("num-samples,n",     prog_opt::value<UINTEGER>()->default_value(100),        "Number of samples for averaging")
        ("tests,t",           prog_opt::value< CLIOptionCommaSeparatedList >()->
                                                  default_value(CLIOptionCommaSeparatedList("all")),               "A comma-separated list of test case ids:" \
                                                                                      "\n0,clc - Cyclic Linear Chain Network" \
                                                                                      "\n1,ca  - Colloidal Dis-/Aggregation Network" \
                                                                                      "\nall - all of the listed above")
        ("sizes,s",           prog_opt::value< CLIOptionCommaSeparatedList >()->
                                                  default_value(CLIOptionCommaSeparatedList("10,100")),            "A comma-separated list of species numbers in the network")
        ("methods,m",         prog_opt::value< CLIOptionCommaSeparatedList >()->
                                                  default_value(CLIOptionCommaSeparatedList("all")),               "A comma-separated list of simulation method ids:" \
                                                                                      "\n0,dm - Gillespie's Direct Method" \
                                                                                      "\n1,pdm - Partial Propensity Direct Method" \
                                                                                      "\n2,pssacr - PSSA with Composition-Rejection Sampling" \
                                                                                      "\n3,spdm - Sorting Partial Propensity Direct Method" \
                                                                                      "\nall - all of the listed above")
        ("verbose,v",                                                                 "Output additional information about the simulation")
        ("quiet,q",                                                                   "Suppress any additional output")
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
   * @copydoc ProgramOptionsBase::resetToDefault()
   */
  virtual void resetToDefault()
  {
    m_bQuiet   = false;
    m_bVerbose = false;

    m_unMethods = 0;
    m_unTests   = 0;

    m_dTimeEnd  = -1.0 * std::numeric_limits<REAL>::infinity();
    m_unSamples = 0;

    m_arunSizes.clear();
  }

  /**
   * @copydoc ProgramOptionsBase::parseVariableMap(prog_opt::variables_map &)
   */
  virtual bool parseVariableMap(prog_opt::variables_map & vm)
  {
    typedef std::map<STRING, UINTEGER> MAPPING_TYPE;
    typedef std::vector<UINTEGER> RESULT_TYPE;

    MAPPING_TYPE mapping;
    RESULT_TYPE  result;

    // output options
    m_bVerbose = (vm.count("verbose") > 0);
    m_bQuiet   = (vm.count("quiet") > 0);

    // simulation parameters
    m_dTimeEnd  = vm["tend"].as<REAL>();
    m_unSamples = vm["num-samples"].as<UINTEGER>();

    m_unTests = 0;
    if(0 == vm.count("tests"))
    {
      m_unTests = tcAll;
    }
    else
    {
      result.clear();
      mapping.clear();

      mapping[STRING("0")] = tcCLC;
      mapping[STRING("clc")] = tcCLC;
      mapping[STRING("1")] = tcCA;
      mapping[STRING("ca")] = tcCA;
      mapping[STRING("all")] = tcAll;

      CLIOptionCommaSeparatedList tests = vm["tests"].as< CLIOptionCommaSeparatedList >();
      tests.parse(mapping, result, true, false, false);

      if(0 == result.size())
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: invalid tests specification. Valid values are:\n\n";
        std::for_each(mapping.begin(), mapping.end(),
                      printPairFirst<MAPPING_TYPE::value_type>(PSSALIB_MPI_CERR_OR_NULL, "\t"));
        PSSALIB_MPI_CERR_OR_NULL << "\n\n";
        return false;
      }
      else
      {
        for(RESULT_TYPE::iterator it = result.begin(); it != result.end(); ++it)
          m_unTests |= *it;
      }
    }

    m_arunSizes.clear();
    if(0 == vm.count("sizes"))
    {
      m_arunSizes.push_back(10);
      m_arunSizes.push_back(100);
    }
    else
    {
      CLIOptionCommaSeparatedList sizes = vm["sizes"].as< CLIOptionCommaSeparatedList >();
      sizes.parse(m_arunSizes, true, false, false);

      if(0 == m_arunSizes.size())
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: network sizes must be positive integers\n\n";
        return false;
      }
    }

    m_unMethods = 0;
    if(0 == vm.count("methods"))
    {
      m_unMethods = pssalib::PSSA::M_All;
    }
    else
    {
      result.clear();
      mapping.clear();

      mapping[STRING("0")] = pssalib::PSSA::M_DM;
      mapping[STRING("dm")] = pssalib::PSSA::M_DM;
      mapping[STRING("1")] = pssalib::PSSA::M_PDM;
      mapping[STRING("pdm")] = pssalib::PSSA::M_PDM;
      mapping[STRING("2")] = pssalib::PSSA::M_PSSACR;
      mapping[STRING("pssacr")] = pssalib::PSSA::M_PSSACR;
      mapping[STRING("3")] = pssalib::PSSA::M_SPDM;
      mapping[STRING("spdm")] = pssalib::PSSA::M_SPDM;
      mapping[STRING("all")] = pssalib::PSSA::M_All;

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
        for(RESULT_TYPE::iterator it = result.begin(); it != result.end(); ++it)
          m_unMethods |= *it;
      }
    }

    return true;
  }

  //////////////////////////////
  // Getters & setters
  //////////////////////////////

  bool isQuietSet() const
  {
    return m_bQuiet;
  }

  bool isVerboseSet() const
  {
    return m_bVerbose;
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

  bool isSamplesSet() const
  {
    return (m_unSamples > 0);
  }

  UINTEGER getSamples() const
  {
    return m_unSamples;
  }

  bool isMethodsSet() const
  {
    return (m_unMethods > 0);
  }

  UINTEGER getMethods() const
  {
    return m_unMethods;
  }

  bool isTestsSet() const
  {
    return (m_unTests > 0);
  }

  UINTEGER getTests() const
  {
    return m_unTests;
  }

  bool isSizesSet() const
  {
    return !m_arunSizes.empty();
  }

  const std::vector< UINTEGER > & getSizes() const
  {
    return m_arunSizes;
  }
};

// Callback for simulation status reporting
void progress_callback(UINTEGER a, UINTEGER b, SHORT c, void * /*user*/)
{
  static SHORT c_old = std::numeric_limits<SHORT>::max();
  if(c != c_old)
  {
    fprintf(stderr, "Progress: sample %u of %u is %hu%% done\n", a, b, c);
    c_old = c;
  }
}

// entry point
int main(int argc, char** argv)
{
  PSSALIB_MPI_IO_INIT;

  try
  {
    // Temporary stream for timing
    STRINGSTREAM ssTiming;
    // Temporary stream for timing
    boost::scoped_array<REAL> arTestCaseTiming;

    UINTEGER methods = 0;

    // Test case object
    Benchmarks benchmarks;
    // Simulation parameters
    pssalib::datamodel::SimulationInfo SimInfo;

    // Parse command line and configuration file, if specified
    prog_opt::variables_map vm;
    if(!benchmarks.processCmdLineArgs(argc, argv, vm)) {
      return -127;
    }

    if(PSSALIB_MPI_IS_MASTER)
    {
      methods = (benchmarks.getMethods() * 0x200040008001ULL & 0x111111111111111ULL) % 0xf;

      arTestCaseTiming.reset(new REAL[2*benchmarks.getSizes().size()*methods]);
      memset(arTestCaseTiming.get(), 0, 2*benchmarks.getSizes().size()*methods*sizeof(REAL));

      if(!ssTiming.good())
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to intialize stream storage for timing info." << std::endl;
        return false;
      }
    }

    ///////////////////////////////////////
    // initialize the SimulationInfo object

    // number of samples
    SimInfo.unSamplesTotal = benchmarks.getSamples();

    // time - 100 time steps
    SimInfo.dTimeEnd  = benchmarks.getTimeEnd();
    SimInfo.dTimeStep = 0.0;

    // output all species
    delete SimInfo.pArSpeciesIds;
    SimInfo.pArSpeciesIds = NULL;

    // suppress all outputs except the desired ones
    SimInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
      | pssalib::datamodel::SimulationInfo::ofTiming
#ifdef PSSALIB_ENGINE_CHECK
      | pssalib::datamodel::SimulationInfo::ofTrajectory
#endif
      ;

    if(benchmarks.isVerboseSet()&&benchmarks.isQuietSet()) {
      PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and 'quiet' flags "
                                  "set, however, the later has priority over the former.\n";
    }

    if(!benchmarks.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofStatus;
    }

    if(benchmarks.isVerboseSet()&&!benchmarks.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog
        | pssalib::datamodel::SimulationInfo::ofInfo
        | pssalib::datamodel::SimulationInfo::ofWarning
        | pssalib::datamodel::SimulationInfo::ofError;
    }

    // redirect streams
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
#else
    SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace
      | pssalib::datamodel::SimulationInfo::eofModuleUpdate
      | pssalib::datamodel::SimulationInfo::eofModuleSampling
      | pssalib::datamodel::SimulationInfo::eofModuleGrouping;

    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cout.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cerr.rdbuf());
#endif
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTiming, ssTiming.rdbuf());

    // create an instance of the simulation engine
    boost::scoped_ptr<pssalib::PSSA> ptrPSSA(new pssalib::PSSA());

    // initialize the call-backs
//     ptrPSSA->SetReactionCallback(&reaction_callback, &gs);
    ptrPSSA->SetProgressCallback(&progress_callback, NULL);

    for(UINTEGER t = 1; t < tcAll; t <<= 1)
    {
      if(0 == (t & benchmarks.getTests()))
        continue;

      if(PSSALIB_MPI_IS_MASTER)
        memset(arTestCaseTiming.get(), 0, 2*benchmarks.getSizes().size()*methods*sizeof(REAL));

      for(UINTEGER k_s = 0; k_s < benchmarks.getSizes().size(); ++k_s)
      {
        // Create an SBML document
        boost::scoped_ptr<LIBSBML_CPP_NAMESPACE::SBMLDocument> pSBMLDoc;
        switch(t)
        {
        case tcCLC:
          pSBMLDoc.reset(CyclicLinearChain::generateSBML(benchmarks.getSizes()[k_s]));
        break;
        case tcCA:
          pSBMLDoc.reset(ColloidalAggregation::generateSBML(benchmarks.getSizes()[k_s]));
        break;
        default:
          PSSALIB_MPI_CERR_OR_NULL << "Error : unknown test case code "
              << t << std::endl;
          return -126;
        break;
        }

        SimInfo.parseSBMLDocument(pSBMLDoc.get());
        pSBMLDoc.reset(NULL);

        for(UINTEGER m = 1, k_m = 0; m < pssalib::PSSA::M_All; m <<= 1)
        {
          if(0 == (m & benchmarks.getMethods()))
            continue;

          // set the simulation method
          if(!ptrPSSA->setMethod((pssalib::PSSA::EMethod)m))
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error : failed to set simulation method to "
              << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << std::endl;
            return -125;
          }

          // run the simulation and collect timing information
          if(ptrPSSA->run(&SimInfo))
          {
            if(PSSALIB_MPI_IS_MASTER)
            {
              // read the timing data from stream
              REAL t; UINTEGER n, k = 0; STRING::value_type c;
              REAL mean = 0.0, M2 = 0.0, delta, delta2, tpr;

              while(ssTiming.good())
              {
                ssTiming >> t >> c >> n;

                ++k;
                tpr = t / REAL(n);

                delta = tpr - mean;
                mean += delta / REAL(k);
                delta2 = tpr - mean;
                M2 += delta * delta2;
              }

              arTestCaseTiming[k_s*methods*2 + 2*k_m] = mean;
              arTestCaseTiming[k_s*methods*2 + 2*k_m + 1] = (M2 > 0.0 ? sqrt(M2 / REAL(k - 1)) : 0.0);

              ssTiming.str(STRING()); ssTiming.clear();
            }
          }
          else
          {
            PSSALIB_MPI_CERR_OR_NULL
              << "FAILED to simulate '" << ptrPSSA->getModelName() << "' using "
              << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << "  ... \n";
            return -124;
          }

          ++k_m;
        }
      }

      // Output results
      if(PSSALIB_MPI_IS_MASTER)
      {
        std::ostream &os = std::cout;
        
        os << "Timing for '" << ptrPSSA->getModelName() << "' over "
           << benchmarks.getSamples() << " samples simulated for "
           << benchmarks.getTimeEnd() << " seconds:\n\n#,";

        for(UINTEGER m = 1; m < pssalib::PSSA::M_All; m <<= 1)
        {
          if(!(m & benchmarks.getMethods()))
            continue;
          os << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << ",";
        }
        os << '\n';

        for(UINTEGER k_s = 0; k_s < benchmarks.getSizes().size(); ++k_s)
        {
          os << benchmarks.getSizes()[k_s] << ",";
          for(UINTEGER k_m = 0; k_m < methods; ++k_m)
          {
            os << arTestCaseTiming[k_s*methods*2 + 2*k_m] << "," << arTestCaseTiming[k_s*methods*2 + 2*k_m + 1] << ",";
          }
          os << '\n';
        }

        os << std::endl;
      }
    }
  }
  catch (prog_opt::error & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: '" << e.what() << "'\n";
    return -5;
  }
  catch (boost::bad_any_cast & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: '" << e.what() << "'\n";
    return -4;
  }
  catch(std::bad_alloc & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << e.what() << ": unable to allocate memory.\n";
    return -3;
  }
  catch(std::runtime_error & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Runtime error : '" << e.what() << "'\n";
    return -2;
  }
  catch(...)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Unexpected exceptional condition encountered.\n";
    return -1;
  }

  return 0;
}
