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
 * Validation of the simulator using analytical results
 */

#include <iostream>

#include "PSSA.h"

#include "util/MPIWrapper.h"
#include "util/ProgramOptionsBase.hpp"

#include "Homoreaction.hpp"
#include "Heteroreaction.hpp"

using namespace pssalib::program_options;

// enumerate all test case ids
enum tagTestCases
{
  tcHomoreaction   = 0x0001,
  tcHeteroreaction = 0x0002,
  tcAll            = 0x0003
};

class Validation : public ProgramOptionsBase
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
  std::vector<UINTEGER> m_arunSamples;

  //! Number of samples
  UINTEGER m_unRepetitions;

//////////////////////////////
// Constructors
public:

  //! Constructor
  Validation()
    : ProgramOptionsBase("Options for Validation")
  {
    resetToDefault();
  }

  //! Destructor
  ~Validation()
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
        ("tests,t",           prog_opt::value< CLIOptionCommaSeparatedList >()->
                                default_value(CLIOptionCommaSeparatedList("all")),    "A comma-separated list of test case ids:"
                                                                                      "\n0,homoreaction - Homoreaction Network"
                                                                                      "\n1,heteroreaction  - Heteroreaction Network"
                                                                                      "\nall - all of the listed above")
        ("repetitions,r",     prog_opt::value<UINTEGER>()->default_value(10),         "Number of repetitions for averaging")
        ("num-samples,n",     prog_opt::value< CLIOptionCommaSeparatedList >()->
                                default_value(CLIOptionCommaSeparatedList(
                                                       "1000,10000,100000")),         "A comma-separated list of sample sizes")
        ("methods,m",         prog_opt::value< CLIOptionCommaSeparatedList >()->
                                default_value(CLIOptionCommaSeparatedList("all")),    "A comma-separated list of simulation method ids:"
                                                                                      "\n0,dm - Gillespie's Direct Method"
                                                                                      "\n1,pdm - Partial Propensity Direct Method"
                                                                                      "\n2,pssacr - PSSA with Composition-Rejection Sampling"
                                                                                      "\n3,spdm - Sorting Partial Propensity Direct Method"
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
    m_unRepetitions = 0;

    m_arunSamples.clear();
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
    m_unRepetitions = vm["repetitions"].as<UINTEGER>();

    if(0 == vm.count("tests"))
    {
      m_unTests = tcAll;
    }
    else
    {
      result.clear();
      mapping.clear();

      mapping[STRING("0")] = tcHomoreaction;
      mapping[STRING("homoreaction")] = tcHomoreaction;
      mapping[STRING("1")] = tcHeteroreaction;
      mapping[STRING("heteroreaction")] = tcHeteroreaction;
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

    m_arunSamples.clear();
    if(0 == vm.count("num-samples"))
    {
      m_arunSamples.push_back(1000);
      m_arunSamples.push_back(10000);
      m_arunSamples.push_back(100000);
    }
    else
    {
      CLIOptionCommaSeparatedList sizes = vm["num-samples"].as< CLIOptionCommaSeparatedList >();
      sizes.parse(m_arunSamples, true, false, false);

      if(0 == m_arunSamples.size())
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: number of samples must be positive integers\n\n";
        return false;
      }
    }

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

  static Hreaction * createTestCase(UINTEGER t)
  {
    switch(t)
    {
    case tcHomoreaction:
      return new Homoreaction;
    break;
    case tcHeteroreaction:
      return new Heteroreaction;
    break;
    default:
      throw std::runtime_error((boost::format("unknown test case id %i") % t).str());
    break;
    }
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
    return !m_arunSamples.empty();
  }

  const std::vector< UINTEGER > & getSamples() const
  {
    return m_arunSamples;
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

  bool isRepetitionsSet() const
  {
    return (m_unRepetitions > 0);
  }

  UINTEGER getRepetitions() const
  {
    return m_unRepetitions;
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
    STRINGSTREAM ssFinalPopulations;
    // Temporary stream for timing
    boost::scoped_array<REAL> arKLDivergence;

    UINTEGER methods = 0;

    // Test case object
    Validation validation;
    // Simulation parameters
    pssalib::datamodel::SimulationInfo SimInfo;

    // Parse command line and configuration file, if specified
    prog_opt::variables_map vm;
    if(!validation.processCmdLineArgs(argc, argv, vm)) {
      return -127;
    }

    if(PSSALIB_MPI_IS_MASTER)
    {
      // count the number of set bits in the lowest byte
      methods = (validation.getMethods() * 0x200040008001ULL & 0x111111111111111ULL) % 0xf;
      arKLDivergence.reset(new REAL[2*validation.getSamples().size()*methods]);
      memset(arKLDivergence.get(), 0, 2*validation.getSamples().size()*methods*sizeof(REAL));

      if(!ssFinalPopulations.good())
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to intialize stream storage for timing info.\n";
        return false;
      }
    }

    ///////////////////////////////////////
    // initialize the SimulationInfo object

    // time - 100 time steps
    SimInfo.dTimeEnd  = validation.getTimeEnd();
    SimInfo.dTimeStep = 0.0;

    // suppress all outputs except the desired ones
    SimInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
#ifndef PSSALIB_ENGINE_CHECK
      | pssalib::datamodel::SimulationInfo::ofFinalPops
#else
      | pssalib::datamodel::SimulationInfo::ofTrajectory
#endif
      ;

    if(validation.isVerboseSet()&&validation.isQuietSet()) {
      PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and "
      "'quiet' flags set, however, the later has priority over the former.\n";
    }

    if(!validation.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofStatus;
    }

    if(validation.isVerboseSet()&&!validation.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog
        | pssalib::datamodel::SimulationInfo::ofInfo
        | pssalib::datamodel::SimulationInfo::ofWarning
        | pssalib::datamodel::SimulationInfo::ofError;
    }

    // redirect streams
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofFinalPops, ssFinalPopulations.rdbuf());
#else
    SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace
      | pssalib::datamodel::SimulationInfo::eofModuleUpdate
      | pssalib::datamodel::SimulationInfo::eofModuleSampling
      | pssalib::datamodel::SimulationInfo::eofModuleGrouping;

    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cout.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cerr.rdbuf());
#endif

    // create an instance of the simulation engine
    boost::scoped_ptr<pssalib::PSSA> ptrPSSA(new pssalib::PSSA());

    // initialize the call-backs
//     ptrPSSA->SetReactionCallback(&reaction_callback, &gs);
    ptrPSSA->SetProgressCallback(&progress_callback, NULL);

    for(UINTEGER t = 1; t < tcAll; t <<= 1)
    {
      if(0 == (t & validation.getTests()))
        continue;

      boost::scoped_ptr<Hreaction> ptrTest(Validation::createTestCase(t));

      // parse the SBML document and initialize the model
      if(!SimInfo.parseSBMLDocument(ptrTest->generateSBML()))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to load the SBML model\n";
        return -126;
      }

      // set species ids for output
      SimInfo.pArSpeciesIds = ptrTest->getSpeciesIds();

      if(PSSALIB_MPI_IS_MASTER)
        memset(arKLDivergence.get(), 0, 2*validation.getSamples().size()*methods*sizeof(REAL));

      for(UINTEGER k_s = 0; k_s < validation.getSamples().size(); ++k_s)
      {
        // number of samples
        const UINTEGER samples = validation.getSamples()[k_s];
        SimInfo.unSamplesTotal = validation.getRepetitions() * samples;

        for(UINTEGER m = 1, k_m = 0; m < pssalib::PSSA::M_All; m <<= 1)
        {
          if(0 == (m & validation.getMethods()))
            continue;

          // set the simulation method
          if(!ptrPSSA->setMethod((pssalib::PSSA::EMethod)m))
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error : failed to set simulation method to "
              << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << std::endl;
            return -126;
          }

          // run the simulation and collect timing information
          if(ptrPSSA->run(&SimInfo))
          {
#ifndef PSSALIB_ENGINE_CHECK
            if(PSSALIB_MPI_IS_MASTER)
            {
              REAL KLdivMean = 0.0, KLdivM2 = 0.0;
              for(UINTEGER r = 0; r < validation.getRepetitions(); ++r)
              {
                // compute the empirical PDF
                std::map<UINTEGER, UINTEGER> mapPDF;
                UINTEGER n, cnt = 0;
                STRING::value_type c;

                while((cnt < samples)&&ssFinalPopulations.good())
                {
                  ssFinalPopulations >> n >> c; cnt++;

                  std::map<UINTEGER, UINTEGER>::iterator
                    itMap = mapPDF.find(n);

                  if(mapPDF.end() == itMap)
                  {
                    mapPDF[n] = 1;
                  }
                  else
                  {
                    itMap->second++;
                  }
                }

                if(cnt < samples)
                {
                  PSSALIB_MPI_CERR_OR_NULL
                    << "FAILED to analyze '" << ptrTest->getName() << "' simulated using "
                    << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m)
                    << ": simulator output contains fewer data points than requested!"
                    << cnt << " < " << samples << std::endl;
                  return -125;
                }

                // output the empirical PDF
                if(0 == r)
                {
                  switch(samples)
                  {
                  case 100:
                  case 1000:
                  case 10000:
                  case 100000:
                  {
                    // output the PDF
                    std::cout << "PDF for '" << ptrTest->getName() << "'  using " << samples << " samples from "
                      << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << std::endl;
                    std::cout << std::setw(4) << "#" << "| Simulated|Analytical|\n";
                    for(std::map<UINTEGER, UINTEGER>::iterator it = mapPDF.begin(); it != mapPDF.end(); it++)
                    {
                      std::cout << std::setw(4) << it->first << ","
                                << std::setw(10) << std::setprecision(5)
                                << ((REAL)it->second)/((REAL)samples) << ","
                                << std::setw(10) << std::setprecision(5)
                                << ptrTest->computePDF(it->first) << std::endl;
                    }
                    std::cout << std::endl;
                    break;
                  }
                  default:
                    break;
                  }
                }

                // compute the divergence
                {
                  const REAL invSamplesSize = 1.0 / REAL(samples),
                            q_eps = 1e-9,    // total probability substitute for all empty bins
                            p_cutoff = 1e-6; // cut-off for the analytical PDF

                  REAL KLdiv = 0.0;
                  UINTEGER i = 0, iBegin = 0, numEmpty = 0, numTotal = 0;
                  while(ptrTest->computePDF(i++) < p_cutoff);
                  iBegin = i;
                  while(true)
                  {
                    REAL p = ptrTest->computePDF(i++);
                    if(p >= p_cutoff)
                    {
                      ++numTotal;

                      std::map<UINTEGER, UINTEGER>::iterator
                        itMap = mapPDF.find(i);

                      if(mapPDF.end() == itMap)
                        ++numEmpty;
                    }
                    else
                      break;
                  }
                  REAL q_corr = REAL(numEmpty) * q_eps / REAL(numTotal - numEmpty);
                  numTotal += iBegin;
                  for(i = iBegin; i < numTotal; ++i)
                  {
                    REAL p = ptrTest->computePDF(i), q;

                    std::map<UINTEGER, UINTEGER>::iterator
                      itMap = mapPDF.find(i);

                    if(mapPDF.end() == itMap)
                      q = q_eps;
                    else
                      q = REAL(itMap->second) * invSamplesSize - q_corr;

                    KLdiv += p * log(p / q);
                  }

                  REAL delta = KLdiv - KLdivMean;
                  KLdivMean += delta / REAL(r + 1);
                  REAL delta2 = KLdiv - KLdivMean;
                  KLdivM2 += delta * delta2;
                }
              }

              arKLDivergence[k_s*methods*2 + 2*k_m] = KLdivMean;
              if((validation.getRepetitions() > 1)&&(KLdivM2 > 0.0))
                arKLDivergence[k_s*methods*2 + 2*k_m + 1] = sqrt(KLdivM2 / REAL(validation.getRepetitions() - 1));

            }

            ssFinalPopulations.str(STRING()); ssFinalPopulations.clear();
#endif
          }
          else
          {
            PSSALIB_MPI_CERR_OR_NULL
              << "FAILED to simulate '" << ptrTest->getName() << "' using "
              << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << "  ... \n";
            return -124;
          }

          ++k_m;
        }
      }

      // Output results
      if(PSSALIB_MPI_IS_MASTER)
      {
#ifndef PSSALIB_ENGINE_CHECK
        std::ostream &os = std::cout;

        os << "Averaged Kullback–Leibler divergence of the simulated PDF from the analytical one for '"
           << ptrTest->getName() << "' over " << validation.getRepetitions() << " repetitions simulated for "
           << validation.getTimeEnd() << " seconds:\n\nSamples,";

        for(UINTEGER m = 1; m < pssalib::PSSA::M_All; m <<= 1)
        {
          if(!(m & validation.getMethods()))
            continue;
          os << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)m) << ",";
        }
        os << '\n';

        for(UINTEGER k_s = 0; k_s < validation.getSamples().size(); ++k_s)
        {
          os << validation.getSamples()[k_s] << ",";
          for(UINTEGER k_m = 0; k_m < methods; ++k_m)
          {
            os << arKLDivergence[k_s*methods*2 + 2*k_m] << "," << arKLDivergence[k_s*methods*2 + 2*k_m + 1] << ",";
          }
          os << '\n';
        }

        os << std::endl;
#endif
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
    PSSALIB_MPI_CERR_OR_NULL << e.what() << ": unable to allocate memory\n";
    return -3;
  }
  catch(std::runtime_error & e)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Runtime error : '" << e.what() << "'\n";
    return -2;
  }
  catch(...)
  {
    PSSALIB_MPI_CERR_OR_NULL << "Unexpected exceptional condition encountered\n";
    return -1;
  }

  return 0;
}
