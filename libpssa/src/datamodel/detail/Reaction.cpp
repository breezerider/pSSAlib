/**
 * @file Reaction.cpp
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
 * Implements a container for certain properties of SBML Reaction objects
 * that are relevant for simulations with Stochastic Simulation Algorithms
 */

#include "../../../include/datamodel/detail/Reaction.h"
#include "../../../include/datamodel/detail/Model.h"
#include "../../../include/datamodel/detail/SpeciesReference.h"
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
  // Copy constructor
  Reaction::Reaction(const Reaction & right)
    : Base(right)
    , m_dFwdRate(right.m_dFwdRate)
    , m_dRevRate(right.m_dRevRate)
    , m_dDelay(right.m_dDelay)
    , m_arSpeciesRefs(NULL)
    , m_unSpeciesRefsCount(right.m_unSpeciesRefsCount)
    , m_unReactants(right.m_unReactants)
    , m_refModel(right.m_refModel)
  {
    if(0 < m_unSpeciesRefsCount)
    {
      m_arSpeciesRefs = new SpeciesReference[m_unSpeciesRefsCount];
      for(UINTEGER i = 0; i < m_unSpeciesRefsCount; ++i)
        m_arSpeciesRefs[i] = right.m_arSpeciesRefs[i];
    }
  }

  /*
   * Release any allocated datastructures
   */
  void Reaction::free()
  {
    if(NULL != m_arSpeciesRefs)
    {
      delete [] m_arSpeciesRefs;
      m_arSpeciesRefs = NULL;
      m_unReactants = m_unSpeciesRefsCount = 0;
    }

    // Call base class method
    Base::free();
  }
#ifdef HAVE_LIBSBML
  /*
   * Assign properties from an SBML Reaction object.
   * 
   * Implementation.
   */
  bool Reaction::assign(const LIBSBML_CPP_NAMESPACE::Reaction *reaction, SBMLHelper & helper)
  {
    // call the base class method
    reinterpret_cast<Base*>(this)->assign(reaction);

    // reversible reaction?
    if(reaction->isSetReversible())
      setReversible(reaction->getReversible());

    // initialize species reference vector
    const LIBSBML_CPP_NAMESPACE::ListOfSpeciesReferences
      *plosrReactants = reaction->getListOfReactants(),
      *plosrProducts = reaction->getListOfProducts();

    std::vector<SpeciesReference> arSpRef; arSpRef.reserve(plosrReactants->size() + plosrProducts->size());
    m_unReactants = 0;

    SpeciesReference spRef(NULL, this);

    if(plosrReactants->size() > 0) 
    {
      for(UINTEGER i = 0; i < plosrReactants->size(); ++i)
      {
        const LIBSBML_CPP_NAMESPACE::SpeciesReference *sr =
          static_cast <const LIBSBML_CPP_NAMESPACE::SpeciesReference *> (plosrReactants->get(i));

        if(!spRef.assign(sr, this))
          return false;

        std::vector<SpeciesReference>::iterator it =
          std::find_if(arSpRef.begin(), arSpRef.end(), spRef);

        if(it != arSpRef.end())
        {
          (*it).setStoichiometry((*it).getStoichiometry() + spRef.getStoichiometry());
        }
        else
        {
          arSpRef.push_back(spRef);
          ++m_unReactants;
        }
      }
    }
    else // reservoir reaction
    {
      spRef.unset(); spRef.makeReservoir();
      arSpRef.push_back(spRef);
      m_unReactants = 1;
    }

    if(plosrProducts->size() > 0) 
    {
      for(UINTEGER i = 0; i < plosrProducts->size(); ++i)
      {
        const LIBSBML_CPP_NAMESPACE::SpeciesReference *sr =
          static_cast <const LIBSBML_CPP_NAMESPACE::SpeciesReference *> (plosrProducts->get(i));

        if(!spRef.assign(sr, this))
          return false;

        std::vector<SpeciesReference>::iterator it =
          std::find_if(arSpRef.begin() + m_unReactants, arSpRef.end(), spRef);

        if(it != arSpRef.end())
        {
          (*it).setStoichiometry((*it).getStoichiometry() + spRef.getStoichiometry());
        }
        else
        {
          arSpRef.push_back(spRef);
        }
      }
    }
    else // reservoir reaction
    {
      spRef.unset(); spRef.makeReservoir();
      arSpRef.push_back(spRef);
    }

    // allocate the species references container & copy the refs
    m_unSpeciesRefsCount = arSpRef.size();
    if(NULL != m_arSpeciesRefs)
      delete [] m_arSpeciesRefs;
    m_arSpeciesRefs = new SpeciesReference[m_unSpeciesRefsCount];
    std::copy(arSpRef.begin(), arSpRef.end(), m_arSpeciesRefs);

    // parse the annotation
    if(!parseAnnotation(reaction, helper))
      return false;

    return true;
  }

  /*
   * Acquire propertie values from the annotation of an SBML Reaction object.
   * 
   * Implementation.
   */
  bool Reaction::parseAnnotation(const LIBSBML_CPP_NAMESPACE::Reaction * reaction, SBMLHelper & helper)
  {
    if(reaction->isSetAnnotation())
    {
      STRING strPrefix("libpSSA");
      XMLNodeTypeDefinition nRoot(STRING(""), STRING(""));
      {
        XMLNodeTypeDefinition nRate(STRING("rate"), strPrefix, XMLCommonTypeDefinition::ctdSingleton | XMLCommonTypeDefinition::ctdMandatory),
                              nFwd(STRING("forward"), strPrefix, XMLCommonTypeDefinition::ctdSingleton | XMLCommonTypeDefinition::ctdMandatory),
                              nRev(STRING("reverse"), strPrefix, this->isReversible() ? (XMLCommonTypeDefinition::ctdSingleton | XMLCommonTypeDefinition::ctdMandatory) : XMLCommonTypeDefinition::ctdSingleton);
        XMLAttributeTypeDefinition aValue(STRING("value"), strPrefix, XMLCommonTypeDefinition::ctdSingleton |XMLCommonTypeDefinition::ctdMandatory);

        nFwd.addAttr(aValue);
        nRev.addAttr(aValue);
        nRate.addChild(nFwd);
        nRate.addChild(nRev);
        nRoot.addChild(nRate);

        XMLNodeTypeDefinition nDelay(STRING("delay"), strPrefix);
        XMLAttributeTypeDefinition aConsuming(STRING("consuming"), strPrefix, XMLCommonTypeDefinition::ctdSingleton),
                                  aNonconsuming(STRING("nonconsuming"), strPrefix, XMLCommonTypeDefinition::ctdSingleton);

        nDelay.addAttr(aValue);
        nDelay.addAttr(aConsuming);
        nDelay.addAttr(aNonconsuming);
        nRoot.addChild(nDelay);
      }

      const LIBSBML_CPP_NAMESPACE::XMLNode * xmlRoot = 
        reaction->getAnnotation();

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

      const XMLNodeTypeDefinition * pnDelay = nRoot.getChild(STRING("delay"), strPrefix),
                                  * pnRate = nRoot.getChild(STRING("rate"), strPrefix),
                                  * pnFwd = (NULL == pnRate) ? NULL : pnRate->getChild(STRING("forward"), strPrefix),
                                  * pnRev = (NULL == pnRate) ? NULL : pnRate->getChild(STRING("reverse"), strPrefix);
      const XMLAttributeTypeDefinition * paValFwd = (NULL == pnFwd) ? NULL : pnFwd->getAttr(STRING("value"), strPrefix),
                                        * paValRev = (NULL == pnRev) ? NULL : pnRev->getAttr(STRING("value"), strPrefix),
                                        * paConsuming = (NULL == pnDelay) ? NULL : pnDelay->getAttr(STRING("consuming"), strPrefix),
                                        * paNonconsuming = (NULL == pnDelay) ? NULL : pnDelay->getAttr(STRING("nonconsuming"), strPrefix),
                                        * paValDelay = (NULL == pnDelay) ? NULL : pnDelay->getAttr(STRING("value"), strPrefix);

      if((NULL == pnDelay)||(NULL == pnRate)||(NULL == pnFwd)||(NULL == pnRev)||(NULL == paValFwd)||(NULL == paValRev)||(NULL == paConsuming)||(NULL == paNonconsuming)||(NULL == paValDelay))
      {
        helper.report(SBMLParserMessage::prtError, 0) << "annotation parser failed: could not initialize parser variables.";
        return false;
      }

      INTEGER nFwdExponent = 1, nRevExponent = 1;
//       UINTEGER unFwdFactor = 1, unRevFactor = 1;
      for(UINTEGER n = 0; n < m_unSpeciesRefsCount; ++n)
      {
        if(n < m_unReactants)
        {
          nFwdExponent -= m_arSpeciesRefs[n].getStoichiometry();
//           unFwdFactor *= gsl_sf_fact(m_arSpeciesRefs[n].getStoichiometryAbs());
          helper.report(SBMLParserMessage::prtTrace, 0) << "reactant " << n << " stoichiometry " << m_arSpeciesRefs[n].getStoichiometry();
        }
        else
        {
          nRevExponent -= m_arSpeciesRefs[n].getStoichiometry();
//           unRevFactor *= gsl_sf_fact(m_arSpeciesRefs[n].getStoichiometryAbs());
          helper.report(SBMLParserMessage::prtTrace, 0) << "product " << (n - m_unReactants) << " stoichiometry " << m_arSpeciesRefs[n].getStoichiometry();
        }
      }

      // Units
      UnitSpecification volumeSpec, substanceSpec, timeSpec;
      // volume specification
      volumeSpec.dExponent = REAL(abs(nFwdExponent) * m_refModel.getCompartmentVolumeDimensions());
      volumeSpec.dMultiplier = REAL(1.0); //pow(m_refModel.getCompartmentVolume(), nFwdExponent);
      // substance specification
      substanceSpec.dExponent = REAL(nFwdExponent);
      substanceSpec.dMultiplier = pow(GSL_CONST_NUM_AVOGADRO, nFwdExponent);
      // time specification
      timeSpec.dExponent = REAL(-1.0);
      timeSpec.dMultiplier = REAL(1.0);

      SBMLHelper::MAP_UNIT_SPEC mapUnits;
      mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE, volumeSpec));
      mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_MOLE, substanceSpec));
      mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_SECOND, timeSpec));

      // Reaction rates
      if(!helper.processValue(paValFwd->getValue(), reaction, m_dFwdRate, mapUnits))
      {
        helper.report(SBMLParserMessage::prtError, paValFwd->getLine()) << "could not process forward rate value.";
        return false;
      }

      if(isReversible())
      {
        // Units
        // volume specification
        volumeSpec.dExponent = REAL(abs(nRevExponent) * m_refModel.getCompartmentVolumeDimensions());
        volumeSpec.dMultiplier = REAL(1.0); //pow(m_refModel.getCompartmentVolume(), nRevExponent);
        // substance specification
        substanceSpec.dExponent = REAL(nFwdExponent);
        substanceSpec.dMultiplier = pow(GSL_CONST_NUM_AVOGADRO, nRevExponent);
        // time specification
        timeSpec.dExponent = REAL(-1.0);
        timeSpec.dMultiplier = REAL(1.0);

        mapUnits.clear();
        mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_METRE, volumeSpec));
        mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_MOLE, substanceSpec));
        mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_SECOND, timeSpec));

        if(!helper.processValue(paValRev->getValue(), reaction, m_dRevRate, mapUnits))
        {
          helper.report(SBMLParserMessage::prtError, paValRev->getLine()) << "could not process reverse rate value.";
          return false;
        }
      }
      else
      {
        if(paValRev->isMatchFound())
        {
          helper.report(SBMLParserMessage::prtWarning, paValRev->getLine()) << "a reverse rate definition was found for an unreversible reaction and was ignored.";
          return false;
        }
      }

      // Reaction delay
      if(pnDelay->isMatchFound())
      {
        if(paValDelay->isMatchFound())
        {
          // Units
          // time specification
          timeSpec.dExponent = REAL(1.0);
          timeSpec.dMultiplier = REAL(1.0);

          mapUnits.clear();
          mapUnits.insert(SBMLHelper::UNIT_SPEC(LIBSBML_CPP_NAMESPACE::UNIT_KIND_SECOND, timeSpec));

          m_unFlags &= rfDelayed;
          if(!helper.processValue(paValDelay->getValue(), reaction, m_dDelay, mapUnits))
          {
            helper.report(SBMLParserMessage::prtError, paValRev->getLine()) << "could not process delay value.";
            return false;
          }

          if(paConsuming->isMatchFound()&&paNonconsuming->isMatchFound())
          {
            helper.report(SBMLParserMessage::prtError, paValDelay->getLine()) << "delay can either be consuming or non-cosuming but both attributes are defined, delay definition ignored.";
            return false;
          }

          if(!paNonconsuming->isMatchFound())
          {
            m_unFlags &= rfConsuming;
          }
        }
        else
        {
          helper.report(SBMLParserMessage::prtWarning, paValDelay->getLine()) << "delay value not defined, delay definition ignored.";
        }
      }
      return true;
    }
    helper.report(SBMLParserMessage::prtError, reaction->getLine()) << "critical error: reaction has no annotation associated with it.";
    return false;
  }
#endif
  /*
    * Allocate species references vector for the reaction.
    * 
    * Implementation
    */
  void Reaction::allocSpeciesRefs(UINTEGER unReactants, UINTEGER unProducts)
  {
    m_unReactants = unReactants;
    m_unSpeciesRefsCount = unReactants + unProducts;
    if(NULL != m_arSpeciesRefs)
      delete [] m_arSpeciesRefs;
    m_arSpeciesRefs = new SpeciesReference[m_unSpeciesRefsCount];
  }

  /*
   * Normalize species references in the reaction.
   * 
   * Implementation
   */
  void Reaction::normalize()
  {
    if(0 == m_unSpeciesRefsCount) return;
    for(UINTEGER i = 0; i < m_unSpeciesRefsCount; ++i)
    {
      SpeciesReference & spI = m_arSpeciesRefs[i];

      bool bReactant = i < m_unReactants;

      for(UINTEGER j = (i + 1); j < (bReactant ? m_unReactants : m_unSpeciesRefsCount); ++j)
      {
        SpeciesReference & spJ = m_arSpeciesRefs[j];

        if(spI(spJ))
        {
          spI.setStoichiometry(spI.getStoichiometry() + spJ.getStoichiometry());
          removeSpeciesReferencesAt(j);
        }
      }
    }
  }
} } } // close namespaces detail, datamodel & pssalib
