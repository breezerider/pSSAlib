/**
 * @file Heteroeaction.hpp
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
 * Implementation of the heteroreaction test case
 */

#include <gsl/gsl_sf_bessel.h>  // Bessel functions
#include <gsl/gsl_sf_gamma.h>   // Factorial
#include <gsl/gsl_sf_exp.h>     // Exponential
#include <gsl/gsl_sf_pow_int.h> // Small integer powers

#include "Hreaction.hpp"

#if __cplusplus > 199711L
#  define CONSTEXPR static constexpr
#else
#  define CONSTEXPR static const
#endif

/**
 * @struct Heteroreaction
 * @brief Simple bimolecular reaction network
 */
struct Heteroreaction : public Hreaction
{
  //! reactor volume
  CONSTEXPR REAL omega = 1.0;
  
  //! rate constants
  CONSTEXPR REAL k1 = 0.04;
  CONSTEXPR REAL k2 = 1.0;

  //! initial population
  CONSTEXPR UINTEGER A0 = 25;
  CONSTEXPR UINTEGER B0 = 1;

  /**
   * @copydoc Hreaction::getName()
   */
  virtual STRING::const_pointer getName() const { return "Heteroreaction"; }

  /**
   * @copydoc Hreaction::generateSBML()
   */
  virtual LIBSBML_CPP_NAMESPACE::SBMLDocument * generateSBML() const
  {
    // Create an SBML document and add a model
    LIBSBML_CPP_NAMESPACE::SBMLDocument * ptrSBMLDocument = new LIBSBML_CPP_NAMESPACE::SBMLDocument(2, 4);

    LIBSBML_CPP_NAMESPACE::Model * ptrSBMLModel =
      ptrSBMLDocument->createModel();

    // set model id
    ptrSBMLModel->setId("Heteroreaction");

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
    compartment->setSize(1.0);
    compartment->setUnits("dimensionless");

    // Species
    LIBSBML_CPP_NAMESPACE::Species * species = NULL;
    // Species A
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("A");
    species->setInitialAmount(A0);
    species->setSubstanceUnits("dimensionless");
    // Species B
    species = ptrSBMLModel->createSpecies();
    species->setCompartment(compartmentName);
    species->setId("B");
    species->setInitialAmount(B0);
    species->setSubstanceUnits("dimensionless");

    // Reactions
    LIBSBML_CPP_NAMESPACE::Reaction * reaction = NULL;
    LIBSBML_CPP_NAMESPACE::SpeciesReference* speciesReference = NULL;
    LIBSBML_CPP_NAMESPACE::Parameter *parameter = NULL;
    LIBSBML_CPP_NAMESPACE::KineticLaw *kineticLaw = NULL;
    boost::format fmtReactionAnnotation("<annotation>\n<libpSSA:rate xmlns:libpSSA=\"uri\">\n"
      "<libpSSA:forward libpSSA:value=\"%s\"/>\n</libpSSA:rate>\n</annotation>");
    //
    // R1: Dimerization reaction
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("A_B_Reaction");
    reaction->setReversible(false);
    // reactant 1
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("A");
    speciesReference->setStoichiometry(1);
    // reactant 2
    speciesReference = reaction->createReactant();
    speciesReference->setSpecies("B");
    speciesReference->setStoichiometry(1);
    // product
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("B");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("kd");
    parameter->setValue(k1 / omega);
    parameter->setUnits(perSecondRate);
    // set reaction annotation
    reaction->setAnnotation((fmtReactionAnnotation % "kd").str());
    //
    // R2: Generation reaction
    reaction = ptrSBMLModel->createReaction();
    reaction->setId("A_Generation");
    reaction->setReversible(false);
    // product
    speciesReference = reaction->createProduct();
    speciesReference->setSpecies("A");
    speciesReference->setStoichiometry(1);
    // rate constant
    kineticLaw = reaction->createKineticLaw();
    parameter = kineticLaw->createParameter();
    parameter->setId("kg");
    parameter->setValue(k2 * omega);
    parameter->setUnits(perSecondRate);
    // set reaction annotation
    reaction->setAnnotation((fmtReactionAnnotation % "kg").str());

    return ptrSBMLDocument;
  }

  /**
   * @copydoc Hreaction::computePDF(UINTEGER)
   */
  virtual REAL computePDF(UINTEGER n) const
  {
    static const REAL K = k2 / k1 * omega * omega / ((REAL)B0);

    return gsl_pow_uint(K, n) / (gsl_sf_exp(K) * gsl_sf_fact(n));
  }
};
