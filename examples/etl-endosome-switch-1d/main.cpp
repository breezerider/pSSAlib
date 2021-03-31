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
 * Early-to-late endosome switch example
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
 * @class ETLSwitch1D
 * @brief Early-to-late endosome switch, 1D model
 */
class ETLSwitch1D : public ProgramOptionsBase
{
//////////////////////////////
// Attributes
private:
  //! units: length - meters, time - seconds

  //! rate constants
  CONSTEXPR REAL k01  = 1.0;
  CONSTEXPR REAL k02  = 1.0e2; // (m^3 / mol)^2 / s
  CONSTEXPR REAL k0m2 = 10.0;
  CONSTEXPR REAL k1   = 1.0;
  CONSTEXPR REAL k21  = 0.1;
  CONSTEXPR REAL k22  = 1.0e2; // (m^3 / mol)^2 / s
  CONSTEXPR REAL k3   = 10.0;

  //! subreactor volume
  CONSTEXPR REAL omega = 4e-21; // m^3

  // initial populations
  REAL S0_t0,    //!< initial population of S0, mol/m^3
       R5_t0,    //!< initial population of R5, mol/m^3
       R7_t0;    //!< initial population of R7, mol/m^3

  REAL D_R5,     //!< diffusion constant for R5
       D_R7;     //!< diffusion constant for R7

  //! Output options
  bool m_bVerbose, m_bQuiet;

  //! Result
  pssalib::datamodel::SimulationInfo::OutputFlags m_eResult;

  //! Method
  pssalib::PSSA::EMethod m_eMethod;

  // simulation end time
  REAL                   m_dTimeEnd;

  //! Number of samples
  UINTEGER m_unSamples;

  // number of subreactors in each direction
  UINTEGER               m_unPoints;

  // output path
  STRING                 m_strOutputPath,
                         m_strFilePattern;

//////////////////////////////
// Constructors
public:

  //! Constructor
  ETLSwitch1D()
    : ProgramOptionsBase("Options for Early-to-late endosome switch (1D) Example")
  {
    resetToDefault();
  }

  //! Destructor
  ~ETLSwitch1D()
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
        ("s0",                prog_opt::value<REAL>()->default_value(0.0),            "initial population of S0")
        ("r5",                prog_opt::value<REAL>()->default_value(0.0),            "initial population of R5")
        ("r7",                prog_opt::value<REAL>()->default_value(0.0),            "initial population of R7")
        ("dr5",               prog_opt::value<REAL>()->default_value(0.0),            "Diffusion rate for species R5")
        ("dr7",               prog_opt::value<REAL>()->default_value(0.0),            "Diffusion rate for species R7")
        ("result,r",          prog_opt::value<STRING>()->required(),                  "Kind of results to output:"
                                                                                      "\n0,panel-c - Produces data points for the phase plot (vary S0(t=0) using --s0 to reproduce the figure)"
                                                                                      "\n0,panel-d - Timing"
                                                                                      "\n0,panel-e - Single trajectory of a reaction-diffusion system")
        ("output-path,o",     prog_opt::value<STRING>()->default_value(
                                                              "output/"),             "Output path")
        ("tend",              prog_opt::value<REAL>()->default_value(0.0),            "End time of the simulation")
        ("num-grid-points,p", prog_opt::value<UINTEGER>()->default_value(0),          "Number of grid points along the stripe")
        ("num-samples,n",     prog_opt::value<UINTEGER>()->default_value(0),          "Number of samples for averaging")
        ("method,m",          prog_opt::value<STRING>()->default_value("spdm"),       "A method id or index:"
                                                                                      "\n0,dm - Gillespie's Direct Method"
                                                                                      "\n1,pdm - Partial Propensity Direct Method"
                                                                                      "\n2,pssacr - PSSA with Composition-Rejection Sampling"
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
      ids.push_back(STRING("R5"));
      ids.push_back(STRING("R7"));
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

    S0_t0= 0.0;
    R5_t0= 0.0;
    R7_t0= 0.0;
    D_R5 = 0.0;
    D_R7 = 0.0;

    m_dTimeEnd = -1.0 * std::numeric_limits<REAL>::infinity();
    m_unSamples = m_unPoints = 0;
    m_strOutputPath = STRING();

    m_eMethod = pssalib::PSSA::M_Invalid;
    m_eResult = pssalib::datamodel::SimulationInfo::ofNone;
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
    S0_t0= vm["s0"].as<REAL>();
    R5_t0= vm["r5"].as<REAL>();
    R7_t0= vm["r7"].as<REAL>();
    D_R5 = vm["dr5"].as<REAL>();
    D_R7 = vm["dr7"].as<REAL>();

    // simulation parameters
    m_dTimeEnd = vm["tend"].as<REAL>();
    m_unSamples = vm["num-samples"].as<UINTEGER>();
    m_unPoints  = vm["num-grid-points"].as<UINTEGER>();
    m_strOutputPath = vm["output-path"].as<STRING>();

    if(0 == vm.count("result"))
    {
      m_eResult = pssalib::datamodel::SimulationInfo::ofNone;
    }
    else
    {
      mapping.clear();

      mapping[STRING("0")] = pssalib::datamodel::SimulationInfo::ofFinalPops;
      mapping[STRING("panel-c")] = pssalib::datamodel::SimulationInfo::ofFinalPops;
      mapping[STRING("1")] = pssalib::datamodel::SimulationInfo::ofTiming;
      mapping[STRING("panel-d")] = pssalib::datamodel::SimulationInfo::ofTiming;
      mapping[STRING("2")] = pssalib::datamodel::SimulationInfo::ofTrajectory;
      mapping[STRING("panel-e")] = pssalib::datamodel::SimulationInfo::ofTrajectory;

      MAPPING_TYPE::iterator it = mapping.find(vm["result"].as<STRING>());

      if(mapping.end() == it)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: invalid result specification.\n";
        return false;
      }
      else
      {
        m_eResult = (pssalib::datamodel::SimulationInfo::OutputFlags)(*it).second;
      }
    }

    switch(m_eResult)
    {
    case pssalib::datamodel::SimulationInfo::ofTrajectory:
    {
      m_unSamples = (m_unSamples > 0 ? m_unSamples : 1);
      m_unPoints = (m_unPoints > 0 ? m_unPoints : 20);
    }
    break;
    case pssalib::datamodel::SimulationInfo::ofTiming:
    {
      m_unSamples = (m_unSamples > 0 ? m_unSamples : 100);
      m_unPoints = (m_unPoints > 0 ? m_unPoints : 1);
    }
    break;
    case pssalib::datamodel::SimulationInfo::ofFinalPops:
    {
      m_unSamples = (m_unSamples > 0 ? m_unSamples : 100);
      m_unPoints = (m_unPoints > 0 ? m_unPoints : 1);
    }
    break;
    default:
      PSSALIB_MPI_CERR_OR_NULL << "Error: unknown result specification\n";
    break;
    }

    if(0 == vm.count("method"))
    {
      m_eMethod = pssalib::PSSA::M_Invalid;
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
        PSSALIB_MPI_CERR_OR_NULL << "Error: invalid method specification.\n";
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
    ptrSBMLModel->setId("ETLSwitch");

    // Unit definitions
    const std::string strRate = "per_second",
                      strVolume = "volume",
                      strSubstance = "substance",
                      strTrimolecular = "per_second_m3_mol_sq";
    LIBSBML_CPP_NAMESPACE_QUALIFIER UnitDefinition * unitDef = NULL;
    LIBSBML_CPP_NAMESPACE_QUALIFIER Unit * unit = NULL;
    // reaction rate
    unitDef =  ptrSBMLModel->createUnitDefinition();
    unitDef->setId(strRate);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_SECOND);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(-1);
    // substance (mole)
    unitDef = ptrSBMLModel->createUnitDefinition();
    unitDef->setId(strSubstance);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_MOLE);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(1.0);
    // volume (m ^ 3)
    unitDef = ptrSBMLModel->createUnitDefinition();
    unitDef->setId(strVolume);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_METRE);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(3.0);
    // rate for trimolecular reactions 1 / s * (m^3 / mole)^2
    unitDef = ptrSBMLModel->createUnitDefinition();
    unitDef->setId(strTrimolecular);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_SECOND);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(-1.0);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_MOLE);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(-2.0);
    unit = unitDef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE_QUALIFIER UNIT_KIND_METRE);
    unit->setMultiplier(1.0); unit->setScale(0.0); unit->setExponent(6.0);

    // Compartment identifier
    const std::string compartmentName = "testVolume";
    // Create a Compartment object
    LIBSBML_CPP_NAMESPACE_QUALIFIER Compartment* compartment =
      ptrSBMLModel->createCompartment();
    compartment->setId(compartmentName);
    compartment->setSize(m_unPoints * omega);
    compartment->setSpatialDimensions((UINTEGER)3);
    compartment->setUnits(strVolume);

    // Species
    LIBSBML_CPP_NAMESPACE_QUALIFIER Species * species = NULL;
    boost::format fmtDiffusionAnnotation("<annotation>\n<libpSSA:diffusion xmlns:libpSSA=\"uri\" libpSSA:value=\"%e\">\n</libpSSA:diffusion>\n</annotation>");
    // Species Rab5
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("R5");
    species->setInitialConcentration(R5_t0);
    species->setSubstanceUnits(strSubstance);
    if((D_R5 > 0.0)&&(LIBSBML_CPP_NAMESPACE_QUALIFIER LIBSBML_OPERATION_SUCCESS != species->setAnnotation((fmtDiffusionAnnotation % D_R5).str())))
      throw std::runtime_error("assigning diffusion annotation to species 'R5' failed!");
    // Species Rab7
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("R7");
    species->setInitialConcentration(R7_t0);
    species->setSubstanceUnits(strSubstance);
    if((D_R7 > 0.0)&&(LIBSBML_CPP_NAMESPACE_QUALIFIER LIBSBML_OPERATION_SUCCESS != species->setAnnotation((fmtDiffusionAnnotation % D_R7).str())))
      throw std::runtime_error("assigning diffusion annotation to species 'R7' failed!");
    // Intermediate Species 0
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("S0");
    species->setInitialConcentration(S0_t0);
    species->setSubstanceUnits(strSubstance);
    // Intermediate Species 1
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("S1");
    species->setInitialConcentration(0);
    species->setSubstanceUnits(strSubstance);

    // Reactions
    LIBSBML_CPP_NAMESPACE_QUALIFIER Reaction * reaction = NULL;
    LIBSBML_CPP_NAMESPACE_QUALIFIER SpeciesReference* speciesReference = NULL;
    LIBSBML_CPP_NAMESPACE_QUALIFIER Parameter *parameter = NULL;
    LIBSBML_CPP_NAMESPACE_QUALIFIER KineticLaw *kineticLaw = NULL;
    boost::format fmtFwdReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"%s\"/>\n</libpSSA:rate>\n</annotation>"),
                                fmtRevReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"%s\"/>\n<libpSSA:reverse libpSSA:value=\"%s\"/>\n</libpSSA:rate>\n</annotation>");
    //
    // R1: Rab5 synthesis
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab5_synthesis");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("R5");
    speciesReference->setStoichiometry(1);
    // product 2
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k01");
    parameter->setValue(k01);
    parameter->setUnits(strRate);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k01").str());
    //
    // R2: Rab5 synthesis - inhibition by Rab7
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab5_synthesis_inhibition");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // reactant 2
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(2);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("S1");
    speciesReference->setStoichiometry(1);
    // product 2
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(2);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k02");
    parameter->setValue(k02);
    parameter->setUnits(strTrimolecular);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k02").str());
    //
    // R3: S1 relaxation
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab5_synthesis_relaxation");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("S1");
    speciesReference->setStoichiometry(1);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k0m2");
    parameter->setValue(k0m2);
    parameter->setUnits(strRate);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k0m2").str());
    //
    // R4: Rab5 degradation
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab5_degradation");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("R5");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k1");
    parameter->setValue(k1);
    parameter->setUnits(strRate);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k1").str());
    //
    // R4: Rab7 synthesis - activation via Rab5
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab7_synthesis_Rab5");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("R5");
    speciesReference->setStoichiometry(1);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(1);
    // product 2
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("R5");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k21");
    parameter->setValue(k21);
    parameter->setUnits(strRate);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k21").str());
    //
    // R5: Rab7 synthesis - inhibition via Rab7
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab7_synthesis_S0");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // reactant 2
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(2);
    // product 1
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("S0");
    speciesReference->setStoichiometry(1);
    // product 2
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(3);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k22");
    parameter->setValue(k22);
    parameter->setUnits(strTrimolecular);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k22").str());
    //
    // R6: Rab7 degradation
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("Rab7_degradation");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("R7");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("k3");
    parameter->setValue(k3);
    parameter->setUnits(strRate);
    // set reaction annotation
    reaction->setAnnotation((fmtFwdReactionAnnotation % "k3").str());

    // pass the document on to the simulation engine
    return SimInfo.parseSBMLDocument(ptrSBMLDocument.release());
  }

  // Population initializer
  static void initialPopulation(pssalib::datamodel::DataModel * ptrData, UINTEGER ** arPtrPopulation, void * etlswitch)
  {
    ETLSwitch1D * ptrETLS = static_cast<ETLSwitch1D *>(etlswitch);
    for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
    {
  //     arPtrPopulation[svi][0] = (svi < ptrETLS->unSubreactors/2 ? ptrETLS->Rab5_t0 : 0.0) * ETLSwitch1D::omega * GSL_CONST_NUM_AVOGADRO;
  //     arPtrPopulation[svi][1] = (svi < ptrETLS->unSubreactors/2 ? 0.0 : ptrETLS->Rab7_t0) * ETLSwitch1D::omega * GSL_CONST_NUM_AVOGADRO;
  //     arPtrPopulation[svi][(svi < ptrETLS->unSubreactors/2 ? 2 : 3)] = ptrETLS->S0_t0 * ETLSwitch1D::omega * GSL_CONST_NUM_AVOGADRO;
      arPtrPopulation[svi][0] = ptrETLS->R5_t0 * ETLSwitch1D::omega * GSL_CONST_NUM_AVOGADRO;
      arPtrPopulation[svi][2] = ptrETLS->S0_t0 * ETLSwitch1D::omega * GSL_CONST_NUM_AVOGADRO;
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
    return (m_unSamples > 0);
  }

  UINTEGER getSamples() const
  {
    return m_unSamples;
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

  bool isResultSet() const
  {
    return (m_eResult != pssalib::datamodel::SimulationInfo::ofNone);
  }

  pssalib::datamodel::SimulationInfo::OutputFlags getResult() const
  {
    return m_eResult;
  }

  REAL getS0() const
  {
    return S0_t0;
  }

  REAL getR5() const
  {
    return R5_t0;
  }

  REAL getR7() const
  {
    return R7_t0;
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
    STRINGSTREAM ssOutput;
    // Test case object
    ETLSwitch1D etlswitch;
    // Simulation parameters
    pssalib::datamodel::SimulationInfo SimInfo;

    // Parse command line and configuration file, if specified
    prog_opt::variables_map vm;
    if(!etlswitch.processCmdLineArgs(argc, argv, vm)) {
      return -127;
    }

    // generate the model
    if(!etlswitch.generateSBML(SimInfo)) {
      return -127;
    }

    ///////////////////////////////////////
    // initialize the SimulationInfo object

    // number of samples
    SimInfo.unSamplesTotal = etlswitch.getSamples();

    // time - 1000 time steps
    SimInfo.dTimeEnd  = etlswitch.getTimeEnd();
    SimInfo.dTimeStep = SimInfo.dTimeEnd / 1e3;

    // output all species
    delete SimInfo.pArSpeciesIds;
    SimInfo.pArSpeciesIds = new std::vector<STRING>(etlswitch.getSpeciesIds());

    // suppress all outputs except the desired ones
    SimInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
#ifndef PSSALIB_ENGINE_CHECK
      | etlswitch.getResult()
#else
      | pssalib::datamodel::SimulationInfo::ofTrajectory
#endif
      ;

    SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace
//       | pssalib::datamodel::SimulationInfo::eofModuleSampling
      | pssalib::datamodel::SimulationInfo::eofModuleGrouping;

    if(etlswitch.isVerboseSet()&&etlswitch.isQuietSet()) {
      PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and 'quiet' flags "
                                  "set, however, the later has priority over the former.\n";
    }

    if(!etlswitch.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofStatus;
    }

    if(etlswitch.isVerboseSet()&&!etlswitch.isQuietSet())
    {
      SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog
        | pssalib::datamodel::SimulationInfo::ofInfo
        | pssalib::datamodel::SimulationInfo::ofWarning
        | pssalib::datamodel::SimulationInfo::ofError;
    }

    // redirect streams
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
    switch(etlswitch.getResult())
    {
    case pssalib::datamodel::SimulationInfo::ofTrajectory:
      SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, ssOutput.rdbuf());
    break;
    case pssalib::datamodel::SimulationInfo::ofTiming:
      if(PSSALIB_MPI_IS_MASTER)
        SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTiming, ssOutput.rdbuf());
    break;
    case pssalib::datamodel::SimulationInfo::ofFinalPops:
      if(PSSALIB_MPI_IS_MASTER)
        SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofFinalPops, ssOutput.rdbuf());
    break;
    default:
      PSSALIB_MPI_CERR_OR_NULL << "Error: unknown result specification\n";
    break;
    }
#else
    SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofTrace
      | pssalib::datamodel::SimulationInfo::eofModuleSampling
      | pssalib::datamodel::SimulationInfo::eofModuleGrouping;

    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cout.rdbuf());
    SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cerr.rdbuf());
#endif

    // set up the domain
#ifndef PSSALIB_ENGINE_CHECK
    SimInfo.setDims(1, (UINTEGER)etlswitch.getPoints());
    SimInfo.eBoundaryConditions = pssalib::datamodel::detail::BC_Periodic;
#endif
    SimInfo.eInitialPopulation = pssalib::datamodel::detail::IP_UserDefined;
    SimInfo.ptrPopulationInitializer = ETLSwitch1D::initialPopulation;
    SimInfo.ptrPopulationInitializerUserData = &etlswitch;

    // create an instance of the simulation engine
    boost::scoped_ptr<pssalib::PSSA> ptrPSSA(new pssalib::PSSA());

    // initialize the call-backs
//     ptrPSSA->SetReactionCallback(&reaction_callback, &etlswitch);
    ptrPSSA->SetProgressCallback(&progress_callback, NULL);

    // set the simulation method
    if(!ptrPSSA->setMethod(etlswitch.getMethod()))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : failed to set simulation method "
        << pssalib::PSSA::getMethodName(etlswitch.getMethod()) << std::endl;
      return -126;
    }

    // run the simulation and collect timing information
    if(ptrPSSA->run(&SimInfo))
    {
#ifndef PSSALIB_ENGINE_CHECK
      if((PSSALIB_MPI_IS_MASTER)||(pssalib::datamodel::SimulationInfo::ofTrajectory == etlswitch.getResult()))
      {
        SimulationDataSource sds;
        if(!sds.load(ssOutput))
        {
          PSSALIB_MPI_CERR_OR_NULL << "Could not load trajectory from the data stream!\n";
          return -125;
        }
        else
        {
          STRING strPath(etlswitch.getOutputPath());
          if(!pssalib::util::makeDir(strPath))
          {
            PSSALIB_MPI_CERR_OR_NULL << "Could not create output path '" << strPath << "'\n";
            return -124;
          }

          // Output
          switch(etlswitch.getResult())
          {
          case pssalib::datamodel::SimulationInfo::ofTrajectory:
          {
            CSVOutputFormatter fmt(STRING(), SimInfo.dTimeStep);

            pssalib::util::makeFilePath(strPath, STRING("subreactor_%i_") + (BOOSTFORMAT("tend=%2.2f_s0=%2.2f_r5=%2.2f_r7=%2.2f.csv") % etlswitch.getTimeEnd() % etlswitch.getS0() % etlswitch.getR5() % etlswitch.getR7()).str(), strPath);
            if(!sds.store(strPath, fmt))
            {
              PSSALIB_MPI_CERR_OR_NULL << "Could not store trajectory as CSV output to '" << strPath << "'\n";
              return -123;
            }
          }
          break;
          case pssalib::datamodel::SimulationInfo::ofTiming:
          {
            CSVOutputFormatter fmt;

            pssalib::util::makeFilePath(strPath, (BOOSTFORMAT("timing_tend=%2.2f_s0=%2.2f_r5=%2.2f_r7=%2.2f_%s.csv") % etlswitch.getTimeEnd() % etlswitch.getS0() % etlswitch.getR5() % etlswitch.getR7() % pssalib::PSSA::getMethodName(etlswitch.getMethod())).str(), strPath);
            if(!sds.store(strPath, fmt))
            {
              PSSALIB_MPI_CERR_OR_NULL << "Could not store timing as CSV output to '" << strPath << "'\n";
              return -123;
            }
          }
          break;
          case pssalib::datamodel::SimulationInfo::ofFinalPops:
          {
            CSVOutputFormatter fmt;

            pssalib::util::makeFilePath(strPath, (BOOSTFORMAT("population_tend=%2.2f_s0=%2.2f_r5=%2.2f_r7=%2.2f.csv") % etlswitch.getTimeEnd() % etlswitch.getS0() % etlswitch.getR5() % etlswitch.getR7()).str(), strPath);
            if(!sds.store(strPath, fmt))
            {
              PSSALIB_MPI_CERR_OR_NULL << "Could not store final population as CSV output to '" << strPath << "'\n";
              return -123;
            }
          }
          break;
          }
        }
      }
#endif
    }
    else
    {
      PSSALIB_MPI_CERR_OR_NULL
        << "FAILED to simulate '" << ptrPSSA->getModelName() << "' using "
        << pssalib::PSSA::getMethodName(etlswitch.getMethod()) << "  ... \n";
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
