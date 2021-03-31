/**
 * @file XMLTypeDefinitions.h
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
 * Defines a container for XML Document Type Definitions (DTD)
 * that are applied to libSBML XMLNode objects during parsing
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_XML_TYPE_DEFINITIONS_H_
#define PSSALIB_DATAMODEL_DETAIL_XML_TYPE_DEFINITIONS_H_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{ 
  // Forward declarations
  class SBMLHelper;

  /**
   * \class XMLCommonTypeDefinition
   * \brief Common Document Type Definitions for XML documents
   * with a predefined structure.
   */
  class XMLCommonTypeDefinition
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagCTDFlags
    {
      ctdSingleton = 0x01, //!< Disallows multiple occurence of this type.
      ctdLastWins = 0x02, //!< If multiple definitions are allowed, gathers the values from last occuring definition, otherwise the first one has priority.
      ctdMandatory = 0x04, //!< If set, this instance is required to find a matching element.

      ctdAllExternal = 0x07,

      ctdMatchFound = 0x08,

      ctdAllInternal = 0x08
    } NodeFlags;

  ////////////////////////////////
  // Attributes
  protected:
    STRING strName,   //!< XML tag name
           strPrefix; //!< XML tag prefix
//            strUri;    //!< XML namespace URI

    //! Line number
    UINTEGER unLineNumber;

    //! Flags
    UINTEGER unFlags;

  /////////////////////////////////////
  // Constructors
  public:
    //! Constructor
    XMLCommonTypeDefinition(const STRING & name, const STRING & prefix, UINTEGER flags = 0)
      : strName(name)
      , strPrefix(prefix)
      , unLineNumber(0)
//       , strUri(uri)
      , unFlags(flags & ctdAllExternal)
    {
      // Do nothing
    };

    //! Destructor
  virtual ~XMLCommonTypeDefinition()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  public:

    /**
     * Get XML tag name.
     * 
     * @return current value
     */
    STRING getName() const
    {
      return strName;
    }

    /**
     * Set the XML tag name.
     * 
     * @param name new value.
     */
    void setName(const STRING & name)
    {
      strName = name;
    }

    /**
     * Get XML tag prefix.
     * 
     * @return current value
     */
    STRING getPrefix() const
    {
      return strPrefix;
    };

    /**
     * Set the XML tag prefix.
     * 
     * @param prefix new value.
     */
    void setPrefix(const STRING & prefix)
    {
      strPrefix = prefix;
    }

    /**
     * Get the number of matching line in the XML file.
     * 
     * @return current value
     */
    UINTEGER getLine() const
    {
      return unLineNumber;
    }

    /**
     * Query if a match has been found.
     * 
     * @return @true if a match was found, false otherwise;
     */
    bool isMatchFound() const
    {
      return (unFlags & ctdMatchFound);
    }

    /**
     * @internal This function is used internally to find a matching instance.
     * @copydoc XMLCommonTypeDefinition::match
     */
    bool operator()(const XMLCommonTypeDefinition & b) const
    {
      return match(b.strPrefix, b.strName);
    }

    /**
     * Matches the prefix and name of this instance to the values
     * provided as arguments.
     * 
     * @return @true if both match exactly, @false otherwise.
     */
    bool match(const STRING & prefix, const STRING & name) const
    {
      return ((strPrefix.compare(prefix)==0)&&(strName.compare(name) == 0));
//       if((strPrefix.compare(prefix)==0)&&(strName.compare(name) == 0))
//       {
//         std::cerr << BOOST_CURRENT_FUNCTION << " : match ('" << prefix
//           << "', '" << name << "')" << std::endl;
//         return true;
//       }
//       else
//       {
//         std::cerr << BOOST_CURRENT_FUNCTION << " : ('" << strPrefix
//           << "', '" << strName << "') DOES NOT match ('" << prefix
//           << "', '" << name << "')" << std::endl;
//         return false;
//       }
    }

    /**
     * Test if current internal state allows to parse a given node.
     * 
     * @return @true if parsing can proceed, @false otherwise.
     */
    bool canParse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, SBMLHelper & helper);

    /**
     * Checks if the parsing did fill out all the required DTD elements.
     * 
     * @return @true if all mandatory elements are filled out, @false otherwise.
     */
  virtual bool finalize(SBMLHelper & helper);
  };

  /**
   * \class XMLAttributeTypeDefinition
   * \brief A DTD-wrapper for XML node attributes.
   * 
   * Represents an attribute of the XML node, carries
   * its traits (i.e., mandatory or not, allowing redefinition).
   */
  class XMLAttributeTypeDefinition : public XMLCommonTypeDefinition
  {
  ////////////////////////////////
  // Attributes
  protected:
    //! Value of the given attribute
    STRING strValue;

  /////////////////////////////////////
  // Constructors
  public:
    //! Constructor
    XMLAttributeTypeDefinition(const STRING & name, const STRING & prefix, UINTEGER flags = 0)//, const STRING & uri = STRING(""))
      : XMLCommonTypeDefinition(name, prefix, flags)//, uri)
      , strValue()
    {
      // Do nothing
    };

    //! Destructor
  virtual ~XMLAttributeTypeDefinition()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  public:

    /**
     * Get XML attribute value.
     * 
     * @return current value
     */
    STRING getValue() const
    {
      return strValue;
    };

    /**
     * Parse attribute value.
     * 
     * @param value Attribute value to parse.
     * @param helper A poiter to XMLAnnotationParser object.
     */
    void parse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, UINTEGER idx, SBMLHelper & helper);
  };

  class XMLNodeTypeDefinition : public XMLCommonTypeDefinition
  {
  ////////////////////////////////
  // Attributes
  protected:

    //! XMLAttributes associated with this node.
    std::vector<XMLAttributeTypeDefinition> colAttributes;

    //! Child nodes.
    std::vector<XMLNodeTypeDefinition> colChildren;

  /////////////////////////////////////
  // Constructors
  public:
    //! Constructor
    XMLNodeTypeDefinition(const STRING & name, const STRING & prefix, UINTEGER flags = 0)
      : XMLCommonTypeDefinition(name, prefix, flags)
    {
      // Do nothing
    };

    //! Destructor
    ~XMLNodeTypeDefinition()
    {
      // Do nothing
    };

    /**************
     * Attributes *
     **************/

    /**
     * Add a new attribute type definition to this node type definition.
     * If another instance with the same name and prefix is present on
     * this node, it will get overwritten by the new one.
     * 
     * @param atd New attribute type definition.
     */
    void addAttr(const XMLAttributeTypeDefinition & atd)
    {
      std::vector<XMLAttributeTypeDefinition>::iterator it =
        std::find_if(colAttributes.begin(), colAttributes.end(), atd);

      if(it != colAttributes.end())
      {
        std::iter_swap(it, colAttributes.rbegin());
        colAttributes.pop_back();
      }
      colAttributes.push_back(atd);
    };

    /**
     * Look up an attribute type definition by its name and prefix.
     * 
     * @param name Name of the sought-for attribute type definition.
     * @param prefix Prefix of the sought-for attribute type definition.
     * @return Pointer to the instance or NULL if none have been found.
     */
    const XMLAttributeTypeDefinition * getAttr(const STRING & name, const STRING & prefix) const
    {
      XMLCommonTypeDefinition t(name, prefix);

      std::vector<XMLAttributeTypeDefinition>::const_iterator it =
        std::find_if(colAttributes.begin(), colAttributes.end(), t);

      if(it != colAttributes.end())
        return &(*it);
      else
        return NULL;
    };

    /**
     * Get an attribute type definition at a given position in the list.
     * 
     * @param index Position of the sought-for attribute type definition.
     * @return Pointer to the instance or NULL if none have been found.
     */
    const XMLAttributeTypeDefinition * getAttr(std::size_t index) const
    {
      std::vector<XMLAttributeTypeDefinition>::const_iterator it =
        colAttributes.begin() + index;

      if(it != colAttributes.end())
        return &(*it);
      else
        return NULL;
    };

    /**
     * Remove an attribute type definition by its name and prefix.
     * 
     * @param name Name of the sought-for attribute type definition.
     * @param prefix Prefix of the sought-for attribute type definition.
     */
    void removeAttr(const STRING & name, const STRING & prefix)
    {
      XMLCommonTypeDefinition t(name, prefix);

      std::vector<XMLAttributeTypeDefinition>::iterator it =
        std::find_if(colAttributes.begin(), colAttributes.end(), t);

      if(it != colAttributes.end())
      {
        std::iter_swap(it, colAttributes.rbegin());
        colAttributes.pop_back();
      }
    };

    /**
     * Remove an attribute type definition at a given position in the list.
     * 
     * @param index Position of the sought-for attribute type definition.
     */
    void removeAttr(std::size_t index)
    {
      std::vector<XMLAttributeTypeDefinition>::iterator it =
        colAttributes.begin() + index;

      if(it != colAttributes.end())
      {
        std::iter_swap(it, colAttributes.rbegin());
        colAttributes.pop_back();
      }
    };

    /************
     * Children *
     ************/

    /**
     * Add a new node type definition to this instance of the class.
     * If another instance with the same name and prefix is present on
     * this node, it will get overwritten by the new one.
     * 
     * @param ntd New node type definition.
     */
    void addChild(const XMLNodeTypeDefinition & ntd)
    {
      std::vector<XMLNodeTypeDefinition>::iterator it =
        std::find_if(colChildren.begin(), colChildren.end(), ntd);

      if(it != colChildren.end())
      {
        std::iter_swap(it, colChildren.rbegin());
        colChildren.pop_back();
      }
      colChildren.push_back(ntd);
    };

    /**
     * Look up a child node type definition by its name and prefix.
     * 
     * @param name Name of the sought-for node type definition.
     * @param prefix Prefix of the sought-for node type definition.
     * @return Pointer to the instance or NULL if none have been found.
     */
    const XMLNodeTypeDefinition * getChild(const STRING & name, const STRING & prefix) const
    {
      XMLCommonTypeDefinition t(name, prefix);

      std::vector<XMLNodeTypeDefinition>::const_iterator it =
        std::find_if(colChildren.begin(), colChildren.end(), t);

      if(it != colChildren.end())
        return &(*it);
      else
        return NULL;
    };

    /**
     * Get a child node type definition at a given position in the list.
     * 
     * @param index Position of the sought-for attribute type definition.
     * @return Pointer to the instance or NULL if none have been found.
     */
    const XMLNodeTypeDefinition * getChild(std::size_t index) const
    {
      std::vector<XMLNodeTypeDefinition>::const_iterator it =
        colChildren.begin() + index;

      if(it != colChildren.end())
        return &(*it);
      else
        return NULL;
    };

    /**
     * Remove a child node type definition by its name and prefix.
     * 
     * @param name Name of the sought-for attribute type definition.
     * @param prefix Prefix of the sought-for attribute type definition.
     */
    void removeChild(const STRING & name, const STRING & prefix)
    {
      XMLCommonTypeDefinition t(name, prefix);

      std::vector<XMLNodeTypeDefinition>::iterator it =
        std::find_if(colChildren.begin(), colChildren.end(), t);

      if(it != colChildren.end())
      {
        std::iter_swap(it, colChildren.rbegin());
        colChildren.pop_back();
      }
    };

    /**
     * Remove a child node type definition at a given position in the list.
     * 
     * @param index Position of the sought-for attribute type definition.
     */
    void removeChild(std::size_t index)
    {
      std::vector<XMLNodeTypeDefinition>::iterator it =
        colChildren.begin() + index;

      if(it != colChildren.end())
      {
        std::iter_swap(it, colChildren.rbegin());
        colChildren.pop_back();
      }
    };

    /**
     * Parse the given XMLNode and all of its siblings according to the 
     * rules defined by the nod and attribute type definitions in this instance.
     * 
     * @param node The XMLNode to parse.
     * @param helper A poiter to XMLAnnotationParser object.
     */
    void parse(const LIBSBML_CPP_NAMESPACE::XMLNode & node, SBMLHelper & helper);

    /**
     * @copydoc XMLCommonTypeDefinition::finalize
     * @override XMLCommonTypeDefinition::finalize
     */
  virtual bool finalize(SBMLHelper & helper);
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_XML_TYPE_DEFEINITIONS_H_ */
