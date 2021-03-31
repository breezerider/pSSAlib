/**
 * @file Model.cpp
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
 * Implements a container for properties of SBML Model object
 * that are relevant for simulating with Stochastic Simulation Algorithms
 */

#include "../../../include/datamodel/detail/Model.h"
#include "../../../include/datamodel/detail/Species.h"
#include "../../../include/datamodel/detail/Reaction.h"
#ifdef HAVE_LIBSBML
#include "../../../include/datamodel/detail/SBMLHelper.hpp"
#endif

#include "../../../include/util/InplaceMemory.h" 

// Implementation
namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /*
   * Release any allocated datastructures
   */
  void Model::free()
  {
    if(m_unFlags & mfShallowCopy)
    {
      m_arSpecies = 0; m_unSpecies = 0;
      m_arReactions = 0; m_unReactions = 0;
      m_unFlags &= ~mfShallowCopy;
    }
    else
    {
      if(NULL != m_arSpecies)
        util::inplace_free(m_unSpecies, &m_arSpecies);

      if(NULL != m_arReactions)
        util::inplace_free(m_unReactions, &m_arReactions);
    }

    // call base class method
    Base::free();
  }

  /*
   * Allocate species vector for the model.
   * 
   * Implementation.
   */
  void Model::allocSpecies(UINTEGER unSpecies)
  {
    if((0 != m_unSpecies)||(NULL != m_arSpecies))
      util::inplace_free(m_unSpecies, &m_arSpecies);
    m_unSpecies = unSpecies;
    m_arSpecies = util::inplace_alloc(m_unSpecies, Species(*this));
  }

  /*
   * Allocate reactions vector for the model.
   * 
   * Implementation.
   */
  void Model::allocReactions(UINTEGER unReactions)
  {
    if((0 != m_unReactions)||(NULL != m_arReactions))
      util::inplace_free(m_unReactions, &m_arReactions);
    m_unReactions = unReactions;
    m_arReactions = util::inplace_alloc(m_unReactions, Reaction(*this));
  }
#ifdef HAVE_LIBSBML
  /*
   * Assign properties from an SBML Model object.
   * 
   * Implementation.
   */
  bool Model::assign(const LIBSBML_CPP_NAMESPACE::Model * model, SBMLHelper & helper, STRING strCompartmentId)
  {
    Base::assign(model);

    const LIBSBML_CPP_NAMESPACE::Compartment * pCompartment;
    if(strCompartmentId.length() > 0)
      pCompartment = model->getCompartment(strCompartmentId);
    else
    {
      pCompartment = model->getCompartment(0);
      strCompartmentId = STRING(pCompartment->getId());
    }
    if(NULL != pCompartment)
    {
      REAL factor = 1.0;
      LIBSBML_CPP_NAMESPACE::UnitDefinition * unitDefSI = NULL;

      if(!pCompartment->isSetUnits())
        helper.report(SBMLParserMessage::prtWarning, pCompartment->getLine()) << "'" << strCompartmentId << "' has no units associated with it.";
      else
      {
        m_uVolumeDims = 0;
        if(NULL != pCompartment->getDerivedUnitDefinition())
        {
          unitDefSI = LIBSBML_CPP_NAMESPACE::UnitDefinition::convertToSI(pCompartment->getDerivedUnitDefinition());

          if(NULL != unitDefSI)
          {
            helper.report(SBMLParserMessage::prtTrace, unitDefSI->getLine()) << "'" << strCompartmentId << "': units '" << LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true) << "'";

            LIBSBML_CPP_NAMESPACE::ListOfUnits * lou = unitDefSI->getListOfUnits();
            for(UINTEGER ui = 0; ui < lou->size(); ++ui)
            {
              LIBSBML_CPP_NAMESPACE::Unit * unit = lou->get(ui);

              if((LIBSBML_CPP_NAMESPACE::UNIT_KIND_DIMENSIONLESS != unit->getKind())&&((LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE != unit->getKind())||(unit->getExponent() < 0.0)))
              {
                helper.report(SBMLParserMessage::prtError, unitDefSI->getLine()) << "'" << strCompartmentId << "': invalid units '" << STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)) << "', can only contain units of length with a positive exponent or be dimensionless!";
                return false;
              }
              else
              {
                if(LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE == unit->getKind())
                  m_uVolumeDims = std::floor(unit->getExponent());
                factor *= unit->getMultiplier() * pow(10, unit->getScale());
              }
            }

            helper.report(SBMLParserMessage::prtTrace, pCompartment->getLine()) << "'" << strCompartmentId << "' has factor = " << factor;
          }
        }
      }

      setCompartmentVolume(pCompartment->getVolume() * factor);

      helper.report(SBMLParserMessage::prtTrace, pCompartment->getLine()) << "'" << strCompartmentId << "' = " << getCompartmentVolume() << " " << ((NULL == unitDefSI) ? " " : STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)));

      if (NULL != unitDefSI)
        delete unitDefSI;
    }
    else
    {
      helper.report(SBMLParserMessage::prtError, 0) << "model has no compartment associated with it.";
      return false;
    }

    // Species: pre-assignment step
    std::vector<const LIBSBML_CPP_NAMESPACE::Species *> arSpeciesPtr;
    const LIBSBML_CPP_NAMESPACE::ListOfSpecies *ploSpecies =
      model->getListOfSpecies();

    for(UINTEGER si = 0; si < ploSpecies->size(); ++si)
    {
      const LIBSBML_CPP_NAMESPACE::Species * ptr = ploSpecies->get(si);
      if(0 == strCompartmentId.compare(STRING(ptr->getCompartment())))
        arSpeciesPtr.push_back(ptr);
    }

    // Reactions: pre-assignment step
    std::vector<const LIBSBML_CPP_NAMESPACE::Reaction *> arReactionsPtr;
    const LIBSBML_CPP_NAMESPACE::ListOfReactions *ploReactions =
      model->getListOfReactions();

    for(UINTEGER ri = 0; ri < ploReactions->size(); ++ri)
    {
      const LIBSBML_CPP_NAMESPACE::Reaction * ptr = ploReactions->get(ri);
      /* FIXME only for SBML v3
       * std::cout << "Reaction #" << ri << " : getCompartment()='" << ptr->getCompartment()
       *   << "'; required='" << strCompartmentId << "'." << std::endl;
       * if(0 == strCompartmentId.compare(STRING(ptr->getCompartment())))
       */
      arReactionsPtr.push_back(ptr);
    }

    // allocate memmory if model contains both reactions and species
    if((0 < arSpeciesPtr.size())&&(0 < arReactionsPtr.size()))
    {
      allocSpecies(arSpeciesPtr.size());
      allocReactions(arReactionsPtr.size());
    }
    else
    {
      free();
      helper.report(SBMLParserMessage::prtError, 0) << "model has no species associated with it.";
      return false;
    }

    // Species: assignment step
    for(UINTEGER si = 0; si < m_unSpecies; ++si)
    {
      if(!assignSpecies(si, arSpeciesPtr[si], helper))
      {
        free();
        helper.report(SBMLParserMessage::prtError, 0) << "species assignment failed.";
        return false;
      }
    }

    // Reactions: assignment step
    for(UINTEGER ri = 0; ri < m_unReactions; ++ri)
    {
      if(!assignReaction(ri, arReactionsPtr[ri], helper))
      {
        free();
        helper.report(SBMLParserMessage::prtError, 0) << "species assignment failed.";
        return false;
      }
    }

    // all OK
    return true;
  }

  /*
   * Add a species to the model.
   * 
   * Implementation.
   */
  bool Model::assignSpecies(UINTEGER i, const LIBSBML_CPP_NAMESPACE::Species * species, SBMLHelper & helper)//, REAL diffusionConstant)
  {
    if(i < m_unSpecies)
    {
      if(!m_arSpecies[i].assign(species, helper))//, this);//, diffusionConstant);
        return false;
      m_mapSpeciesId2Index.insert( MAP_ID2IDX::value_type(species->getId(), i) );
      return true;
    }
    else
      return false;
  }

  /*
   * Add a reaction to the model.
   * 
   * Implementation.
   */
  bool Model::assignReaction(UINTEGER i, const LIBSBML_CPP_NAMESPACE::Reaction * reaction, SBMLHelper & helper)//, REAL fwdRate, REAL revRate, REAL delay, bool consuming)
  {
    if(i < m_unReactions)
      return m_arReactions[i].assign(reaction, helper);
    else
      return false;
  }
#endif
  /*
   * Get a species from the model.
   * 
   * Implementation.
   */
  Species * Model::getSpecies(UINTEGER i) const
  {
    if(i < m_unSpecies)
      return &m_arSpecies[i];
    else
      return NULL;
  }

  /*
   * Get a species from the model.
   * 
   * Implementation.
   */
  UINTEGER Model::getSpeciesIndex(const Species * species) const
  {
    if((species > m_arSpecies)&&(species < (m_arSpecies + m_unSpecies)))
      return species - m_arSpecies;
    else
    {
      STRINGSTREAM ssTemp;
      ssTemp << BOOST_CURRENT_FUNCTION << " : invalid species pointer '" << species << "'";
      throw std::runtime_error(ssTemp.str());
    }
  }

  /*
   * Get a species from the model.
   * 
   * Implementation
   */
  Reaction * Model::getReaction(UINTEGER i) const
  {
    if(i < m_unReactions)
      return &m_arReactions[i];
    else
      return NULL;
  }

  /*
   * Normalize reactions in the model.
   */
  void Model::normalize()
  {
    if(NULL == m_arReactions) return;
    // Reactions: assignment step
    for(UINTEGER i = 0; i < m_unReactions; ++i)
        m_arReactions[i].normalize();
  }

} } } // close namespaces detail, datamodel & pssalib
