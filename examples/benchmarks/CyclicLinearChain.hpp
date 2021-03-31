/**
 * @file TestCyclicLinearChain.cpp
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
 * Implementation of the Cyclic Linear Chain test case.
 */

/**
 * @struct CyclicLinearChain
 * @brief Idealized weakly coupled reaction network
 */
struct CyclicLinearChain
{
  /**
   * Generate the SBML model
   * 
   * @return an @c SBMLDocument object representing the network
   */
  static LIBSBML_CPP_NAMESPACE::SBMLDocument * generateSBML(unsigned int unSpeciesNum)
  {
    // Create an SBML document and add a model
    LIBSBML_CPP_NAMESPACE::SBMLDocument * ptrSBMLDocument = new LIBSBML_CPP_NAMESPACE::SBMLDocument(2, 4);

    LIBSBML_CPP_NAMESPACE::Model * ptrSBMLModel =
      ptrSBMLDocument->createModel();

    ptrSBMLModel->setId("CyclicLinearChain");

    // Unit definitions
    const std::string perSecondRate = "per_second";
    // A UnitDefinition for reaction rates
    LIBSBML_CPP_NAMESPACE::UnitDefinition *unitdef =
      ptrSBMLModel->createUnitDefinition();
    unitdef->setId(perSecondRate);
    // CThe respective unit object
    LIBSBML_CPP_NAMESPACE::Unit *unit =
      unitdef->createUnit();
    unit->setKind(LIBSBML_CPP_NAMESPACE::UNIT_KIND_SECOND);
    unit->setExponent(-1);

    // Compartment identifier
    const std::string compartmentName = "testVolume";
    // Create a Compartment object
    LIBSBML_CPP_NAMESPACE::Compartment* compartment =
      ptrSBMLModel->createCompartment();
    compartment->setId(compartmentName);
    compartment->setSize(1);
    compartment->setUnits("dimensionless");

    // Species
    boost::format fmtSpeciesName("S%d");
    for(UINTEGER s = 0; s < unSpeciesNum; s++)
    {
      LIBSBML_CPP_NAMESPACE::Species * species =
        ptrSBMLModel->createSpecies();
      species->setCompartment(compartmentName);
      species->setId((fmtSpeciesName % s).str());
      species->setInitialAmount(1.0);
      species->setSubstanceUnits("dimensionless");
    }

    // Reactions
    boost::format fmtReactionName("R%d");
    boost::format fmtReactionRate("k%d");
    boost::format fmtReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"k%d\"/>\n</libpSSA:rate>\n</annotation>");
    LIBSBML_CPP_NAMESPACE::Reaction * reaction = NULL;
    LIBSBML_CPP_NAMESPACE::SpeciesReference* speciesReference = NULL;
    LIBSBML_CPP_NAMESPACE::Parameter *parameter = NULL;
    LIBSBML_CPP_NAMESPACE::KineticLaw *kineticLaw = NULL;
    for(UINTEGER s = 0; s < unSpeciesNum - 1; s++)
    {
      reaction = ptrSBMLModel->createReaction();
      reaction->setId((fmtReactionName % s).str());
      reaction->setReversible(false);

      // reactant
      speciesReference = reaction->createReactant();
      speciesReference->setSpecies((fmtSpeciesName % s).str());

      // product
      speciesReference = reaction->createProduct();
      speciesReference->setSpecies((fmtSpeciesName % (s + 1)).str());

      // rate constant
      kineticLaw = reaction->createKineticLaw();
      parameter = kineticLaw->createParameter();
      parameter->setId((fmtReactionRate % s).str());
      parameter->setValue(1.0);
      parameter->setUnits(perSecondRate);

      // set reaction annotation
      reaction->setAnnotation((fmtReactionAnnotation % s).str());
    }

    // close the loop
    reaction = ptrSBMLModel->createReaction();
    reaction->setId((fmtReactionName % (unSpeciesNum - 1)).str());
    reaction->setReversible(false);

    // reactant
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies((fmtSpeciesName % (unSpeciesNum - 1)).str());

    // product
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies((fmtSpeciesName % 0).str());

    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId((fmtReactionRate % (unSpeciesNum - 1)).str());
    parameter->setValue(1.0);
    parameter->setUnits(perSecondRate);

    // set reaction annotation
    reaction->setAnnotation((fmtReactionAnnotation % (unSpeciesNum - 1)).str());

    return ptrSBMLDocument;
  }
};
