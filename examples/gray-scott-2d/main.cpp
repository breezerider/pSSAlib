/**
 * @file main.cpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @version 1.0.0
 * @date Mon, 10 Feb 2017
 * @section LICENSE
 * 
 * The GNU LGPL v3 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Gray-Scott 2D Example
 */

#include <iostream>

#include <boost/program_options.hpp>

#include "PSSA.h"

#include "util/MPIWrapper.h"
#include "util/FileSystem.h"
#include "util/ProgramOptionsBase.hpp"
#include "util/SimulationDataSource.hpp"

using namespace pssalib::program_options;

#if __cplusplus > 199711L
#  define CONSTEXPR static constexpr
#else
#  define CONSTEXPR static const
#endif

/**
 * @class GrayScott2D
 * @brief Gary-Scott 2D system
 */
class GrayScott2D : public ProgramOptionsBase
{
//////////////////////////////
// Attributes
private:

  CONSTEXPR REAL H = 0.01;

  //! Output options
  bool m_bVerbose, m_bQuiet;

  // model parameters
  REAL F,
       k,
       k1,
       u;

  REAL D_A, //!< diffusion constant for A
       D_B; //!< diffusion constant for B

  //! Method
  pssalib::PSSA::EMethod m_eMethod;

  // simulation end time
  REAL                   m_dTimeEnd;
  
  // number of subreactors in each direction
  UINTEGER               m_unPoints;

  // output path
  STRING                 m_strOutputPath,
                         m_strFilePattern;

//////////////////////////////
// Constructors
public:

  //! Constructor
  GrayScott2D()
    : ProgramOptionsBase("Options for Gray-Scott 2D Example")
  {
    resetToDefault();
  }

  //! Destructor
  ~GrayScott2D()
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
        ("F",                 prog_opt::value<REAL>()->default_value(0.043),          "Parameter F")
        ("k",                 prog_opt::value<REAL>()->default_value(0.065),          "Parameter k")
        ("k1",                prog_opt::value<REAL>()->default_value(1.0),            "Parameter k1")
        ("u",                 prog_opt::value<REAL>()->default_value(1e7),            "Parameter u")
        ("da",                prog_opt::value<REAL>()->default_value(8e9),            "Diffusion rate for species A")
        ("db",                prog_opt::value<REAL>()->default_value(4e9),            "Diffusion rate for species B")
        ("output-path,o",     prog_opt::value<STRING>()->default_value(
                                                         "vtk/vanilla/"),             "Output path")
        ("file-pattern,f",    prog_opt::value<STRING>()->default_value(
                                                      "sequence_%i.vtk"),             "Output file pattern (use %i for frame number)")
        ("tend",              prog_opt::value<REAL>()->default_value(1.0),            "End time of the simulation")
        ("num-grid-points,p", prog_opt::value<UINTEGER>()->default_value(64),         "Number of grid points in each direction of a square lattice")
        ("method,m",          prog_opt::value<STRING>(),                              "A method id or index:" \
                                                                                      "\n0,dm - Gillespie's Direct Method" \
                                                                                      "\n1,pdm - Partial Propensity Direct Method" \
                                                                                      "\n2,pssacr - PSSA with Composition-Rejection Sampling" \
                                                                                      "\n3,spdm - Sorting Partial Propensity Direct Method")
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

  // Get ids of species in the current model
  const std::vector<STRING> & getSpeciesIds() const
  {
    static std::vector<STRING> ids;
    if(ids.empty())
    {
      ids.push_back(STRING("S0"));
      ids.push_back(STRING("S1"));
    }
    return ids;
  }

  /**
   * @copydoc ProgramOptionsBase::resetToDefault()
   */
  virtual void resetToDefault()
  {
    m_bQuiet   = false;
    m_bVerbose = false;

    F   = 0.043;
    k   = 0.069;
    k1  = 1.0;
    u   = 1e6;
    D_A = 8e7;
    D_B = 4e7;

    m_dTimeEnd = -1.0 * std::numeric_limits<REAL>::infinity();
    m_unPoints = 0;
    m_strFilePattern = m_strOutputPath = STRING();

    m_eMethod = pssalib::PSSA::M_Invalid;
  }

  /**
   * @copydoc ProgramOptionsBase::parseVariableMap(prog_opt::variables_map &)
   */
  virtual bool parseVariableMap(prog_opt::variables_map & vm)
  {
    typedef std::map<STRING, UINTEGER> MAPPING_TYPE;

    MAPPING_TYPE mapping;

    // output options
    m_bVerbose = (vm.count("verbose") > 0);
    m_bQuiet   = (vm.count("quiet") > 0);

    // model parameters
    F = vm["F"].as<REAL>();
    k = vm["k"].as<REAL>();
    k1 = vm["k1"].as<REAL>();
    u = vm["u"].as<REAL>();
    D_A = vm["da"].as<REAL>();
    D_B = vm["db"].as<REAL>();

    // simulation parameters
    m_dTimeEnd = vm["tend"].as<REAL>() * 2e3/u/u;
    m_unPoints  = vm["num-grid-points"].as<UINTEGER>();
    m_strOutputPath = vm["output-path"].as<STRING>();
    m_strFilePattern = vm["file-pattern"].as<STRING>();

    if(0 == vm.count("method"))
    {
      m_eMethod = pssalib::PSSA::M_PDM;
    }
    else
    {
      mapping.clear();

      mapping[STRING("0")] = pssalib::PSSA::M_DM;
      mapping[STRING("dm")] = pssalib::PSSA::M_DM;
      mapping[STRING("1")] = pssalib::PSSA::M_PDM;
      mapping[STRING("pdm")] = pssalib::PSSA::M_PDM;
      mapping[STRING("2")] = pssalib::PSSA::M_PSSACR;
      mapping[STRING("pssacr")] = pssalib::PSSA::M_PSSACR;
      mapping[STRING("3")] = pssalib::PSSA::M_SPDM;
      mapping[STRING("spdm")] = pssalib::PSSA::M_SPDM;
      
      MAPPING_TYPE::iterator it = mapping.find(vm["method"].as<STRING>());
      
      if(mapping.end() == it)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: invalid method specification.";
        return false;
      }
      else
      {
        m_eMethod = (pssalib::PSSA::EMethod)(*it).second;
      }
    }

    return true;
  }
  
  /**
   * Generates an SBML model using current parameter values.
   * @return @true if parser succeeds, @false otherwise.
   */
  bool generateSBML(pssalib::datamodel::SimulationInfo & SimInfo)
  {
    // Create an SBML document and add a model
    std::auto_ptr<LIBSBML_CPP_NAMESPACE::SBMLDocument> ptrSBMLDocument
      (new LIBSBML_CPP_NAMESPACE::SBMLDocument(2, 4));

    LIBSBML_CPP_NAMESPACE::Model * ptrSBMLModel =
      ptrSBMLDocument->createModel();

    // id
    ptrSBMLModel->setId("GrayScott");

    // Compartment identifier
    const std::string compartmentName = "testVolume";
    // Create a Compartment object
    LIBSBML_CPP_NAMESPACE::Compartment* compartment =
      ptrSBMLModel->createCompartment();
    compartment->setId(compartmentName);
    compartment->setSize(REAL(m_unPoints * m_unPoints) * GrayScott2D::H * GrayScott2D::H);
    compartment->setSpatialDimensions((UINTEGER)2);
//     compartment->setUnits(strVolume);

    // Species
    LIBSBML_CPP_NAMESPACE::Species * species = NULL;
    boost::format fmtDiffusionAnnotation("<annotation>\n<libpSSA:diffusion xmlns:libpSSA=\"uri\" "
      "libpSSA:value=\"%e\">\n</libpSSA:diffusion>\n</annotation>");
    //
    // Species A
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("A");
    species->setInitialAmount(0);
    // set species annotation
    species->setAnnotation((fmtDiffusionAnnotation % (D_A)).str());
    //
    // Species B
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("B");
    species->setInitialAmount(0);
    // set species annotation
    species->setAnnotation((fmtDiffusionAnnotation % (D_B)).str());

    // Reactions
    LIBSBML_CPP_NAMESPACE::Reaction * reaction = NULL;
    LIBSBML_CPP_NAMESPACE::SpeciesReference* speciesReference = NULL;
    boost::format fmtFwdReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"%e\"/>\n</libpSSA:rate>\n</annotation>");
    //
    // R1: A synthesis
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("A_synthesis");
    reaction->setReversible(false);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("A");
    speciesReference->setStoichiometry(1);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % (F*k1*u*u*u)).str());
    //
    // R2: A degradation
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("A_degradation");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("A");
    speciesReference->setStoichiometry(1);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % (F*k1*u*u)).str());
    //
    // R3: A to B conversion
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("A_conversion");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("B");
    speciesReference->setStoichiometry(2);
    // reactant 2
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("A");
    speciesReference->setStoichiometry(1);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("B");
    speciesReference->setStoichiometry(3);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % (k1)).str());
    //
    // R4: B degradation
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("B_degradation");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("B");
    speciesReference->setStoichiometry(1);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % ((k + F)*k1*u*u)).str());

    // pass the document on to the simulation engine
    return SimInfo.parseSBMLDocument(ptrSBMLDocument.release());
  }

  /**
   * Population initializer (see @file typedef.h for argument definition).
   */
  static void initialPopulation(pssalib::datamodel::DataModel * ptrData, UINTEGER ** arPtrPopulation, void * grayscott)
  {
    GrayScott2D * ptrGS = static_cast<GrayScott2D *>(grayscott);
    gsl_rng * ptrRNG = gsl_rng_alloc(gsl_rng_default);

    const REAL UHH = ptrGS->u * GrayScott2D::H * GrayScott2D::H;
    const UINTEGER NN = ptrGS->getPoints();
    const UINTEGER lo = std::floor(0.375*NN);
    const UINTEGER hi = std::floor(0.625*NN);

    for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
    {
      REAL r = (REAL)gsl_rng_uniform (ptrRNG);

      UINTEGER a = svi % NN;
      UINTEGER b = svi / NN;
#ifndef PSSALIB_ENGINE_CHECK
      if (a > lo && a < hi && b > lo && b < hi)
      {
        arPtrPopulation[svi][0] = UHH/2.0 + (0.04*(r-0.5)*UHH + 0.5);
        arPtrPopulation[svi][1] = UHH/4.0 + (0.02*(r-0.5)*UHH + 0.5);
      }
      else
      {
        arPtrPopulation[svi][0] = UHH;
        arPtrPopulation[svi][1] = 0;
      }
#else
      arPtrPopulation[svi][0] = UHH;
      arPtrPopulation[svi][1] = UHH;
#endif
    }

    gsl_rng_free(ptrRNG);
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

  bool isPointsSet() const
  {
    return (m_unPoints > 0);
  }

  UINTEGER getPoints() const
  {
    return m_unPoints;
  }

  bool isOutputPathSet() const
  {
    return !m_strOutputPath.empty();
  }

  STRING getOutputPath() const
  {
    return m_strOutputPath;
  }

  bool isFilePatternSet() const
  {
    return !m_strFilePattern.empty();
  }

  STRING getFilePattern() const
  {
    return m_strFilePattern;
  }

  bool isMethodSet() const
  {
    return (m_eMethod != pssalib::PSSA::M_Invalid);
  }

  pssalib::PSSA::EMethod getMethod() const
  {
    return m_eMethod;
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
    // Temporary stream for trajectory
    STRINGSTREAM ssTrajectory;
    // Test case object
    GrayScott2D grayscott;
    // Simulation parameters
    pssalib::datamodel::SimulationInfo SimInfo;

    // Parse command line arguments and configuration file (if specified)
    prog_opt::variables_map vm;
    if(!grayscott.processCmdLineArgs(argc, argv, vm)) {
      return -127;
    }

    // generate the model
    if(!grayscott.generateSBML(SimInfo)) {
      return -127;
    }

    ///////////////////////////////////////
    // initialize the SimulationInfo object

    // number of samples
    SimInfo.unSamplesTotal = 1;

    // time - 100 time steps
    SimInfo.dTimeEnd  = grayscott.getTimeEnd();
    SimInfo.dTimeStep = SimInfo.dTimeEnd / 100.0;

    // output all species
    delete SimInfo.pArSpeciesIds;
    SimInfo.pArSpeciesIds = NULL;

    // suppress all outputs except the desired ones
    SimInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
      | pssalib::datamodel::SimulationInfo::ofTrajectory;

    if(grayscott.isVerboseSet()&&grayscott.isQuietSet()) {
      PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and 'quiet' flags "
                                  "set, however, the later has priority over the former.\n";
    }

    if(!grayscott.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofStatus;
    }

    if(grayscott.isVerboseSet()&&!grayscott.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog
        | pssalib::datamodel::SimulationInfo::ofInfo
        | pssalib::datamodel::SimulationInfo::ofWarning
        | pssalib::datamodel::SimulationInfo::ofError;
    }

    // redirect streams
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, ssTrajectory.rdbuf());
#else
    SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace
      | pssalib::datamodel::SimulationInfo::eofModuleSampling
      | pssalib::datamodel::SimulationInfo::eofModuleGrouping;

    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cout.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cerr.rdbuf());
#endif

    // set up the domain
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setDims(2, grayscott.getPoints(), grayscott.getPoints());
    SimInfo.eBoundaryConditions = pssalib::datamodel::detail::BC_Periodic;
#endif
    SimInfo.eInitialPopulation = pssalib::datamodel::detail::IP_UserDefined;
    SimInfo.ptrPopulationInitializer = GrayScott2D::initialPopulation;
    SimInfo.ptrPopulationInitializerUserData = &grayscott;

    // create an instance of the simulation engine
    boost::scoped_ptr<pssalib::PSSA> ptrPSSA(new pssalib::PSSA());

    // initialize the call-backs
//     ptrPSSA->SetReactionCallback(&reaction_callback, &grayscott);
    ptrPSSA->SetProgressCallback(&progress_callback, NULL);

    // set the simulation method
    if(!ptrPSSA->setMethod(grayscott.getMethod()))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to set simulation method "
        << pssalib::PSSA::getMethodName(grayscott.getMethod()) << std::endl;
      return -126;
    }

    // run the simulation and collect timing information
    if(ptrPSSA->run(&SimInfo))
    {
#ifndef PSSALIB_ENGINE_CHECK
      SimulationDataSource sds;
      if(!sds.load(ssTrajectory))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Could not load trajectory from the data stream!\n";
        return -125;
      }
      else
      {
        VTKOutputFormatter fmt(SimInfo.getDimsCount(), SimInfo.getDims(), grayscott.getSpeciesIds());
        STRING strPath(grayscott.getOutputPath());
        if(!pssalib::util::makeDir(strPath))
        {
          PSSALIB_MPI_CERR_OR_NULL << "Could not create output path '" << strPath << "'\n";
          return -124;
        }
        pssalib::util::makeFilePath(strPath, grayscott.getFilePattern(), strPath);
        if(!sds.store(strPath, fmt))
        {
          PSSALIB_MPI_CERR_OR_NULL << "Could not store trajectory as VTK output to '" << strPath << "'\n";
          return -123;
        }
      }
#endif
    }
    else
    {
      PSSALIB_MPI_CERR_OR_NULL
        << "FAILED to simulate '" << ptrPSSA->getModelName() << "' using "
        << pssalib::PSSA::getMethodName(grayscott.getMethod()) << "  ... \n";
      return -122;
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
