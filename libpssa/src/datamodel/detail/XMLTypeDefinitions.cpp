/**
 * @file XMLTypeDefinitions.cpp
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
 * Implements a container for XML Document Type Definitions (DTD)
 * that are applied to libSBML XMLNode objects during parsing
 */

#ifdef HAVE_CONFIG_H
  #include "config.h" // current build configuration
#endif

#ifdef HAVE_LIBSBML

#include "../../../include/datamodel/detail/XMLTypeDefinitions.h"
#include "../../../include/datamodel/detail/SBMLHelper.hpp"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /*
   * Test if current internal state allows to parse a given node.
   * 
   * Implementation.
   */
  bool XMLCommonTypeDefinition::canParse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, SBMLHelper & helper)
  {
    // respect single definition requirement
    if(unFlags & ctdMatchFound)
    {
      if(unFlags & ctdSingleton)
      {
        helper.report(SBMLParserMessage::prtWarning, node.getLine()) << " multiple definitions of node '" << node.getPrefix() << ":" << node.getName() << "' are ignored.";
        return false;
      }

      if(unFlags & ctdLastWins)
      {
        // TODO unset all values
      }
    }

    // mark match
    unFlags |= ctdMatchFound;

    // store line number
    unLineNumber = node.getLine();

    return true;
  }

  /*
   * Checks if the parsing did fill out all the required DTD elements.
   * 
   * Implementation.
   */
  bool XMLCommonTypeDefinition::finalize(SBMLHelper & helper)
  {
    if((unFlags & ctdMandatory)&&(0 == (unFlags & ctdMatchFound)))
    {
      helper.report(SBMLParserMessage::prtError, unLineNumber) << "entity '" << strPrefix << ":" << strName << "' required by DTD not found.";
      return false;
    }

    return true;
  }

  /**
   * Parse attribute value.
   * 
   * Implementation.
   */
  void XMLAttributeTypeDefinition::parse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, UINTEGER idx, SBMLHelper & helper)
  {
    if(!canParse(node, helper))
      return;

    helper.report(SBMLParserMessage::prtTrace, node.getLine()) << "node '" << node.getPrefix() << ":" << node.getName() << "' has attribute '" <<node.getAttrPrefix(idx) << ":" << node.getAttrName(idx) << " with value '" << node.getAttrValue(idx) << "' assigned.";

    strValue = node.getAttrValue(idx);
  }

  /*
   * Parse the given XMLNode and all of its siblings according to the 
   * rules defined by the nod and attribute type definitions in this instance.
   * 
   * Implementation.
   */
  void XMLNodeTypeDefinition::parse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, SBMLHelper & helper)
  {
    if(!canParse(node, helper))
      return;

    helper.report(SBMLParserMessage::prtTrace, node.getLine()) << "node '" << node.getPrefix() << ":" << node.getName() << "' is parsed by '" << strPrefix << ":" << strName << "' next.";

    if(0 != colAttributes.size())
    {
      if(0 == node.getAttributesLength())
      {
        helper.report(SBMLParserMessage::prtWarning, node.getLine()) << "no attributes are found for node '" << node.getPrefix() << ":" << node.getName() << "', however, attributes are present in the document type definition.";
      }
      else
      {
        for(INTEGER ai = 0; ai < node.getAttributesLength(); ++ai)
        {
          XMLCommonTypeDefinition t(node.getAttrName(ai), node.getAttrPrefix(ai));

          std::vector<XMLAttributeTypeDefinition>::iterator it =
            std::find_if(colAttributes.begin(), colAttributes.end(), t);

          if(it != colAttributes.end())
            (*it).parse(node, ai, helper);
          else
            helper.report(SBMLParserMessage::prtWarning, node.getLine()) << "unrecognized attribute '" << node.getAttrPrefix(ai) << ":" << node.getAttrName(ai) << "' in node '" << node.getPrefix() << ":" << node.getName() << "'.";
        }
      }
    }
    else
    {
      if(0 != node.getAttributesLength())
      {
        helper.report(SBMLParserMessage::prtWarning, node.getLine()) << "attributes are found for node '" << node.getPrefix() << ":" << node.getName() << "', however, none are required by the document type definition.";
      }
    }

    if(0 != colChildren.size())
    {
      if(0 == node.getNumChildren())
      {
        helper.report(SBMLParserMessage::prtWarning, node.getLine()) << "no children are found for node '" << node.getPrefix() << ":" << node.getName() << "', however, child nodes are present in the document type definition.";
      }
      else
      {
        for(UINTEGER ci = 0; ci < node.getNumChildren(); ++ci)
        {
          const LIBSBML_CPP_NAMESPACE::XMLNode & cn = node.getChild(ci);
          XMLCommonTypeDefinition t(cn.getName(), cn.getPrefix());

          std::vector<XMLNodeTypeDefinition>::iterator it =
            std::find_if(colChildren.begin(), colChildren.end(), t);

          if(it != colChildren.end())
            (*it).parse(cn, helper);
          else
            helper.report(SBMLParserMessage::prtWarning, cn.getLine()) << "unrecognized child node '" << cn.getPrefix() << ":" << cn.getName() << "' in node '" << node.getPrefix() << ":" << node.getName() << "'.";
        }
      }
    }
    else
    {
      if(0 != node.getNumChildren())
      {
        helper.report(SBMLParserMessage::prtWarning, node.getLine()) << "children nodes are found for node '" << node.getPrefix() << ":" << node.getName() << "', however, none are required by the document type definition.";
      }
    }
  }

  /*
   * Implementation.
   */
  bool XMLNodeTypeDefinition::finalize(SBMLHelper & helper)
  {
    bool resultSelf = XMLCommonTypeDefinition::finalize(helper);

    if(!resultSelf)
      return false;
    else if(0 == (unFlags & ctdMatchFound))
      return true; // ignore absence of non-mandatory nodes

    bool resultAttr = true;
    if(0 != colAttributes.size())
    {
      for(std::vector<XMLAttributeTypeDefinition>::iterator it = colAttributes.begin();
          it != colAttributes.end(); ++it)
      {
        resultAttr = resultAttr && (*it).finalize(helper);
      }
    }

    bool resultChild = true;
    if(0 != colChildren.size())
    {
      for(std::vector<XMLNodeTypeDefinition>::iterator it = colChildren.begin();
          it != colChildren.end(); ++it)
      {
        resultChild = resultChild && (*it).finalize(helper);
      }
    }

    return resultSelf && resultAttr && resultChild;
  }

} } } // close namespaces detail, datamodel & pssalib

#else

#warning "SBML support has been disabled!"

#endif
