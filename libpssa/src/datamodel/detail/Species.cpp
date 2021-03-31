/**
 * @file Species.cpp
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
 * Implements a container for certain properties of SBML Species objects
 * that are relevant for simulations with Stochastic Simulation Algorithms
 */

#include "../../../include/datamodel/detail/Species.h"
#include "../../../include/datamodel/detail/Model.h"
#ifdef HAVE_LIBSBML
#include "../../../include/datamodel/detail/XMLTypeDefinitions.h"
#include "../../../include/datamodel/detail/SBMLHelper.hpp"
#endif

namespace pssalib
{
namespace datamodel
{
namespace detail
{
#ifdef HAVE_LIBSBML
  /*
   * Assign properties from an SBML Species object.
   * 
   * Implementation
   */
  bool Species::assign(const LIBSBML_CPP_NAMESPACE::Species *species, SBMLHelper & helper)
  {
    // call the base class method
    reinterpret_cast<Base*>(this)->assign(species);

    // parse annotation
    if(!parseAnnotation(species, helper))
      return false;

    // get units
    REAL factor = 1.0;
    LIBSBML_CPP_NAMESPACE::UnitDefinition * unitDefSI = LIBSBML_CPP_NAMESPACE::UnitDefinition::convertToSI(species->getDerivedUnitDefinition());

    if(species->isSetInitialConcentration())
    {
      if(NULL != unitDefSI)
      {
        LIBSBML_CPP_NAMESPACE::ListOfUnits * lou = unitDefSI->getListOfUnits();
        for(UINTEGER ui = 0; ui < lou->size(); ++ui)
        {
          LIBSBML_CPP_NAMESPACE::Unit * unit = lou->get(ui);
          switch(unit->getKind())
          {
          case LIBSBML_CPP_NAMESPACE::UNIT_KIND_MOLE:
            factor *= pow(GSL_CONST_NUM_AVOGADRO, unit->getExponent());
          break;
          case LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE:
          case LIBSBML_CPP_NAMESPACE::UNIT_KIND_DIMENSIONLESS:
            // Do nothing
          break;
          default:
          {
            helper.report(SBMLParserMessage::prtError, unitDefSI->getLine()) << "Units for species initial amounts/concentrations may only contain 'mole', 'metre' or 'dimensionless'. However, '" << STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)) << "' found for '" << this->toString() << "'.";
            return false;
          }
          break;
          }

          factor *= pow(unit->getMultiplier() * pow(10, unit->getScale()), unit->getExponent());
        }
      }
      else
      {
        helper.report(SBMLParserMessage::prtWarning, species->getLine()) << "'" << this->toString() << "': cannot obtain units.";
      }
    }

    // initial amount
    REAL amount = 0.0;
    if(species->isSetInitialAmount())
      amount = species->getInitialAmount();
    else if(species->isSetInitialConcentration())
    {
      factor *= m_refModel.getCompartmentVolume();
      amount = species->getInitialConcentration();
    }
    else
    {
      helper.report(SBMLParserMessage::prtWarning, species->getLine()) << "'" << this->toString() << "': both initial concentration and amount are not set. Assuming none present.";
    }

    helper.report(SBMLParserMessage::prtTrace, species->getLine()) << "'" << this->toString() << "': initial amount = " << amount << " " << ((!species->isSetInitialConcentration())||(NULL == unitDefSI) ? " " : STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true))) << " ==> " << round(amount * factor);

    if (NULL != unitDefSI)
      delete unitDefSI;

    setInitialAmount(round(amount * factor));

    // boundary conditions
    setConstant(species->isSetConstant()&&species->getConstant());
    setBoundaryCondition(species->isSetBoundaryCondition()&&species->getBoundaryCondition());

    return true;
  }

  /*
   * Acquire propertie values from the annotation of an SBML Species object.
   * 
   * Implementation.
   */
  bool Species::parseAnnotation(const LIBSBML_CPP_NAMESPACE::Species * species, SBMLHelper & helper)
  {
    if(species->isSetAnnotation())
    {
      STRING strPrefix("libpSSA");
      XMLNodeTypeDefinition nRoot(STRING(""), STRING(""));
      {
        XMLNodeTypeDefinition nDiffusion(STRING("diffusion"), strPrefix, XMLCommonTypeDefinition::ctdSingleton);
        XMLAttributeTypeDefinition aValue(STRING("value"), strPrefix, XMLCommonTypeDefinition::ctdSingleton);

        nDiffusion.addAttr(aValue);

        nRoot.addChild(nDiffusion);
      }

      const LIBSBML_CPP_NAMESPACE::XMLNode * xmlRoot = 
        species->getAnnotation();

      // this should never happen
      if(NULL == xmlRoot)
      {
        helper.report(SBMLParserMessage::prtError, 0) << "could not obtain the annotation object from the SBML object.";
        return false;
      }

      nRoot.parse(*xmlRoot, helper);
      if(!nRoot.finalize(helper))
      {
        helper.report(SBMLParserMessage::prtError, 0) << "annotation parser failed: could not initialize a helper.";
        return false;
      }

      const XMLNodeTypeDefinition * pnDiffusion = nRoot.getChild(STRING("diffusion"), strPrefix);
      const XMLAttributeTypeDefinition * paValDiffusion = (NULL == pnDiffusion) ? NULL : pnDiffusion->getAttr(STRING("value"), strPrefix);

      if((NULL == pnDiffusion)||(NULL == paValDiffusion))
      {
        helper.report(SBMLParserMessage::prtError, 0) << "annotation parser failed: could not initialize parser variables.";
        return false;
      }

      // Species diffusion
      if(pnDiffusion->isMatchFound())
      {
        if(paValDiffusion->isMatchFound())
        {
          // Units
          UnitSpecification volumeSpec, substanceSpec, timeSpec;
          // volume specification
          volumeSpec.dExponent = REAL(2.0);
          // time specification
          timeSpec.dExponent = REAL(-1.0);

          SBMLHelper::MAP_UNIT_SPEC mapUnits;
          mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE, volumeSpec));
          mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_SECOND, timeSpec));

          // Parse diffusion constant value
          REAL tempDiff = 0.0;
          if(!helper.processValue(paValDiffusion->getValue(), species, tempDiff, mapUnits))
          {
            helper.report(SBMLParserMessage::prtError, paValDiffusion->getLine()) << "could not process diffusion rate value.";
            return false;
          }

          helper.report(SBMLParserMessage::prtTrace, paValDiffusion->getLine()) << "diffusion rate value = " << tempDiff;

          // Assign it to the species
          setDiffusionConstant(tempDiff);
        }
        else
        {
          helper.report(SBMLParserMessage::prtWarning, paValDiffusion->getLine()) << "diffusion value not defined, diffusion definition ignored.";
        }
      }

      return true;
    }
    else
      return true; // no annotation is acceptable
    return false;
  }
#endif
  /*
   * Get the index of this species in the model.
   * 
   * Implementation
   */
  UINTEGER Species::getIndex() const
  {
    return m_refModel.getSpeciesIndex(this);
  }

} } } // close namespaces detail, datamodel & pssalib
