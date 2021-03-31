/**
 * @file ColloidalAggregation.hpp
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
 * Implementation of the Colloidal Aggregation test case.
 */

/**
 * @struct ColloidalAggregation
 * @brief Idealized strongly coupled reaction network
 */
struct ColloidalAggregation
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

    // set model id
    ptrSBMLModel->setId("ColloidalAggregation");

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

    boost::format fmtSpeciesName("S%d");
    for(UINTEGER s = 0; s < unSpeciesNum; s++)
    {
      LIBSBML_CPP_NAMESPACE::Species * species =
        ptrSBMLModel->createSpecies();
      species->setCompartment(compartmentName);
      species->setId((fmtSpeciesName % s).str());
      species->setInitialAmount(1.0);
  //       if(0 == s)
  //         species->setInitialAmount(unSpeciesNum);
  //       else
  //         species->setInitialAmount(0);
      species->setSubstanceUnits("dimensionless");
    }

    boost::format fmtReactionName("R%d");
    boost::format fmtReactionRate("k%d");
    boost::format fmtReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"%1.1f\"/>\n</libpSSA:rate>\n</annotation>");
    LIBSBML_CPP_NAMESPACE::Reaction * reaction = NULL;
    LIBSBML_CPP_NAMESPACE::SpeciesReference* speciesReference = NULL;
    LIBSBML_CPP_NAMESPACE::Parameter *parameter = NULL;
    LIBSBML_CPP_NAMESPACE::KineticLaw *kineticLaw = NULL;

    UINTEGER r = 0;

    // S_n + S_m --> S_(n+m)
    for(UINTEGER s1 = 0; s1 < (unSpeciesNum / 2); s1++)
    {
      for(UINTEGER s2 = s1; s2 < (unSpeciesNum - s1 - 1); s2++, r++)
      {
        reaction = ptrSBMLModel->createReaction();
        reaction->setId((fmtReactionName % r).str());
        reaction->setReversible(false);

        // reactant 1
        speciesReference = reaction->createReactant();
        speciesReference->setSpecies((fmtSpeciesName % s1).str());

        // reactant 2
        speciesReference = reaction->createReactant();
        speciesReference->setSpecies((fmtSpeciesName % s2).str());

        // product
        speciesReference = reaction->createProduct();
        speciesReference->setSpecies((fmtSpeciesName % (s1 + s2 + 1)).str());

  //         // rate constant
  //         kineticLaw = reaction->createKineticLaw();
  //         parameter = kineticLaw->createParameter();
  //         parameter->setId((fmtReactionRate % r).str());
  //         parameter->setValue(1.0);
  //         parameter->setUnits(perSecondRate);

        // set reaction annotation
        reaction->setAnnotation((fmtReactionAnnotation % 1.0).str());
      }
    }

    // S_p --> S_q + S_(p-q)
    for(UINTEGER s1 = 1; s1 < unSpeciesNum; s1++)
    {
      for(UINTEGER s2 = 0; s2 < ((s1 - 1) / 2 + 1); s2++, r++)
      {
        reaction = ptrSBMLModel->createReaction();
        reaction->setId((fmtReactionName % r).str());
        reaction->setReversible(false);

        // reactant
        speciesReference = reaction->createReactant();
        speciesReference->setSpecies((fmtSpeciesName % s1).str());

        // product 1
        speciesReference = reaction->createProduct();
        speciesReference->setSpecies((fmtSpeciesName % s2).str());

        // product 2
        speciesReference = reaction->createProduct();
        speciesReference->setSpecies((fmtSpeciesName % (s1 - s2 - 1)).str());

  //         // rate constant
  //         kineticLaw = reaction->createKineticLaw();
  //         parameter = kineticLaw->createParameter();
  //         parameter->setId((fmtReactionRate % r).str());
  //         parameter->setValue(1.0);
  //         parameter->setUnits(perSecondRate);

        // set reaction annotation
        reaction->setAnnotation((fmtReactionAnnotation % 1.0).str());
      }
    }

    return ptrSBMLDocument;
  }
};
