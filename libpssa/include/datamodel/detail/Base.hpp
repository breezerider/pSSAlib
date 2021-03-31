/**
 * @file Base.hpp
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
 * Declares a datatype containing for chemical species
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_BASE_HPP_
#define PSSALIB_DATAMODEL_DETAIL_BASE_HPP_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /**
   * @class Base
   * @brief Declares a container for basic SBML node properties.
   */
  class Base
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagBaseFlags
    {
      bfIdSet = 0x01,
      bfNameSet = bfIdSet << 1,

      bfAll = bfIdSet + bfNameSet
    } BaseFlags;

  ////////////////////////////////
  // Attributes
  protected:
    // Basic properties
    STRING m_strId,     //!< 'id' property
           m_strName;   //!< 'name' property

    UINTEGER m_unFlags; //!<  Flags

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    Base()
      : m_strId()
      , m_strName()
      , m_unFlags(0)
    {
      // Do nothing
    };
#ifdef HAVE_LIBSBML
    //! Constructor
  explicit Base(const LIBSBML_CPP_NAMESPACE::SBase * base)
      : m_strId()
      , m_strName()
      , m_unFlags(0)
    {
      assign(base);
    };
#endif
    //! Copy constructor
    Base(const Base & right)
      : m_strId(right.m_strId)
      , m_strName(right.m_strName)
      , m_unFlags(right.m_unFlags)
    {
      // Do nothing
    };

    //! Destructor
  virtual ~Base()
    {
      free();
    };

  /////////////////////////////////////
  // Methods
  public:
    /**
     * Release any allocated datastructures
     */
  virtual void free()
    {
      // Do nothing
    };

    /**
     * Reset all properties' values.
     */
  virtual void unset()
    {
      free();

      m_strId.clear();
      m_strName.clear();
      m_unFlags = 0;
    };
#ifdef HAVE_LIBSBML
    /**
     * Assign properties from an SBML NamedSBase object.
     * 
     * @param base an SBML SBase object.
     */
    void assign(const LIBSBML_CPP_NAMESPACE::SBase * base)
    {
      unset(); // clear properties

      if(NULL != base)
      {
        // id
        if(base->isSetId())
          setId(base->getId());

        // name
        if(base->isSetName())
          setName(base->getName());
      }
    };
#endif
    /**
     * Get the SBML id.
     * 
     * @return current value.
     */
  inline const STRING getId() const
    {
      if(m_unFlags & bfIdSet)
        return m_strId;
      return STRING();
    };

    /**
     * Set the SBML id.
     * 
     * @param id new value.
     */
  inline void setId(const STRING & id)
    {
      m_unFlags |= bfIdSet;
      m_strId = id;
    };

    /**
     * Get the SBML name.
     * 
     * @return current value.
     */
  inline STRING getName() const
    {
      if(m_unFlags & bfNameSet)
        return m_strName;
      return STRING();
    };

    /**
     * Set the SBML name.
     * 
     * @param name new value.
     */
  inline void setName(const STRING & name)
    {
      m_unFlags |= bfNameSet;
      m_strName = name;
    };

    /**
     * Get a string represantation of this object.
     * 
     * @return string representing this object.
     */
    const STRING & toString() const
    {
      static STRING result;
      if(m_unFlags & bfNameSet)
        result = m_strName + " [" + m_strId + "]";
      else
        result = "[" + m_strId + "]";
      return result;
    }

    /**
     * Copy attribute values from another instance
     */
  virtual void copy(const Base & other)
    {
      m_unFlags = other.m_unFlags;
      m_strId = other.m_strId;
      m_strName = other.m_strName;
    }

    /**
     * Swap members with another instance
     */
  virtual void swap(Base & other)
    {
      std::swap(m_unFlags, other.m_unFlags);
      std::swap(m_strId, other.m_strId);
      std::swap(m_strName, other.m_strName);
    }
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_BASE_HPP_ */
