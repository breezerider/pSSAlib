/**
 * @file Reaction.h
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
 * Declares a container for properties of SBML Reaction objects that are
 * relevant for simulations using Stochastic Simulation Algorithms
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_REACTION_H_
#define PSSALIB_DATAMODEL_DETAIL_REACTION_H_

#include "Base.hpp"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  // Forward declarations
  class Model;
  class SpeciesReference;
  class SBMLHelper;

  /**
   * @class Reaction
   * @brief Declares a datatype for chemical reactions
   */
  class Reaction : public Base
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagReactionFlags
    {
      rfReversible = (bfAll + 1),
      rfForward = rfReversible << 1,

      rfDelayed = rfReversible << 2,
      rfConsuming = rfReversible << 3,

      rfAll = rfReversible + rfForward + rfDelayed + rfConsuming
    } ReactionFlags;

  ////////////////////////////////
  // Attributes
  protected:
    // Rate constants
    REAL m_dFwdRate, //!< Forward rate
         m_dRevRate; //!< Reverse rate

    //! Time delay for delayed reactions
    REAL m_dDelay;

    //! Reactants vector
    SpeciesReference *m_arSpeciesRefs;

    UINTEGER m_unSpeciesRefsCount,//!< Total number of species refs
             m_unReactants;       //!< Number of reactants

    //! Container
    Model const & m_refModel;

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    Reaction(Model const & model)
      : Base()
      , m_dFwdRate(std::numeric_limits<REAL>::min())
      , m_dRevRate(std::numeric_limits<REAL>::min())
      , m_dDelay(std::numeric_limits<REAL>::min())
      , m_arSpeciesRefs(NULL)
      , m_unSpeciesRefsCount(0)
      , m_unReactants(0)
      , m_refModel(model)
    {
      // Do nothing
    };

    //! Copy constructor
    Reaction(const Reaction & right);

    //! Destructor
    ~Reaction()
    {
      free();
    };

  /////////////////////////////////////
  // Methods
  private:
#ifdef HAVE_LIBSBML
    /**
     * Acquire propertie values from the annotation of an SBML Reaction object.
     * 
     * @param reaction an SBML @link Reaction object.
     */
    bool parseAnnotation(const LIBSBML_CPP_NAMESPACE::Reaction * reaction, SBMLHelper & helper);
#endif
  public:

    /**
     * Release any allocated datastructures
     */
  virtual void free();

    /**
     * Reset all properties' values.
     */
  virtual void unset()
    {
      // call base class method
      Base::unset();

      m_dFwdRate = std::numeric_limits<REAL>::min();
      m_dRevRate = std::numeric_limits<REAL>::min();
      m_dDelay = std::numeric_limits<REAL>::min();
    };
#ifdef HAVE_LIBSBML
    /**
     * Assign properties from an SBML Reaction object.
     * 
     * @param reaction an SBML @link Reaction object.
     * @param helper an @link SBMLHelper object.
     * @return @true if assignment was successful, @false otherwise.
     */
    bool assign(const LIBSBML_CPP_NAMESPACE::Reaction * reaction, SBMLHelper & helper);
#endif

    /**
     * Allocate species references vector for the reaction.
     * 
     * @param unReactants number of reactants.
     * @param unProducts number of products.
     */
    void allocSpeciesRefs(UINTEGER unReactants, UINTEGER unProducts);

    /**
     * Get the Model container object for this instance.
     * 
     * @return Model object containing this instance.
     */
    const Model * getModel() const { return &m_refModel; };

    /**
     * Query whether this reaction is reversible.
     * 
     * @return @true if this reaction is reversible, @false otherwise.
     */
  inline void setReversible(bool reversible) { reversible ? (m_unFlags |= rfReversible) : (m_unFlags &= ~rfReversible); };

    /**
     * Query whether this reaction is reversible.
     * 
     * @return @true if this reaction is reversible, @false otherwise.
     */
  inline bool isReversible() const { return (m_unFlags & rfReversible); };

    /**
     * Get the forward reaction rate.
     * 
     * @return current value.
     */
  inline REAL getForwardRate() const { return m_dFwdRate; };

    /**
     * Set the forward reaction rate.
     * 
     * @param dC new value.
     */
  inline void setForwardRate(const REAL dC) { m_dFwdRate = dC; };

    /**
     * Get the reverse reaction rate.
     * 
     * @return current value.
     */
  inline REAL getReverseRate() const
    {
      if(m_unFlags & rfReversible)
        return m_dRevRate;
      return std::numeric_limits<REAL>::min();
    };

    /**
     * Set the reverse reaction rate.
     * 
     * @param dC new value.
     */
  inline void setReverseRate(const REAL dC)
    {
      if(m_unFlags & rfReversible)
        m_dRevRate = dC;
    };

    /**
     * Get the number of reactant species references associated with this reaction.
     * 
     * @return current value.
     */
 inline UINTEGER getReactantsCount() const { return m_unReactants; };

    /**
     * Get the number of product references associated with this reaction.
     * 
     * @return current value.
     */
  inline UINTEGER getProductsCount() const { return m_unSpeciesRefsCount - m_unReactants; };

    /**
     * Get the total number of species references associated with this reaction.
     * 
     * @return current value.
     */
  inline UINTEGER getSpeciesReferencesCount() const { return m_unSpeciesRefsCount; };

    /**
     * Get a reactant species reference from the respective list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the reactants
     * list or NULL if the index exceeds the upper bound.
     */
  inline const SpeciesReference * getReactantsListAt(UINTEGER n) const;

    /**
     * Get a reactant species reference from the respective list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the reactants
     * list or NULL if the index exceeds the upper bound.
     */
  inline SpeciesReference * getReactantsListAt(UINTEGER n);

    /**
     * Get a product species reference from the respective list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the products
     * list or NULL if the index exceeds the upper bound.
     */
  inline const SpeciesReference * getProductsListAt(UINTEGER n) const;

    /**
     * Get a product species reference from the respective list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the products
     * list or NULL if the index exceeds the upper bound.
     */
  inline SpeciesReference * getProductsListAt(UINTEGER n);

    /**
     * Get a species reference from the list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the
     * list or NULL if the index exceeds the upper bound.
     */
  inline const SpeciesReference * getSpeciesReferenceAt(UINTEGER n) const;

    /**
     * Get a species reference from the list.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the
     * list or NULL if the index exceeds the upper bound.
     */
  inline SpeciesReference * getSpeciesReferenceAt(UINTEGER n);

    /**
     * Get the species reference list.
     * 
     * @return Array of @link getSpeciesReferencesCount() SpeciesReference objects.
     */
  inline const SpeciesReference * getSpeciesReferences() const;

    /**
     * Remove the species reference at the given positions in the list.
     * 
     * @param n Index of the species reference to remove.
     * @return @true if parameters are within list bounds, @false otherwise.
     */
  inline bool removeSpeciesReferencesAt(UINTEGER n);

    /**
     * Swap the species reference at the given positions in the list.
     * 
     * @param n1 Index of the first species reference to exchange places with the other.
     * @param n2 Index of the second species reference to exchange places with the other.
     * @return @true if parameters are within list bounds, @false otherwise.
     */
  inline bool swapSpeciesReferencesAt(UINTEGER n1, UINTEGER n2);

  /**
   * Normalize species references in the reaction.
   */
  void normalize();

    /**
     * Get the reaction delay.
     * 
     * @return current value.
     */
  inline UINTEGER getDelay() const
    {
      if(m_unFlags & rfDelayed)
        return m_dDelay;
      return 0.0;
    };

    /**
     * Set the reaction delay. A negative value resets it to 
     * an ordinary reaction.
     * 
     * @param dD new value.
     */
  inline void setDelay(const REAL dD)
    {
      if(dD > 0.0)
      {
        m_unFlags |= rfDelayed;
        m_dDelay = dD;
      }
      else
      {
        m_unFlags &= ~rfDelayed;
        m_dDelay = std::numeric_limits<REAL>::min();
      }
    };

    /**
     * Check if the reaction delay is set.
     * 
     * @return @true if property is set, @false otherwise.
     */
  inline bool isSetDelay() const { return (bool)(m_unFlags & rfDelayed); };

    /**
     * Set this delay type to consuming.
     * 
     * @param consuming new value.
     */
  inline void setDelayConsuming(bool consuming) { consuming ? (m_unFlags |= rfConsuming) : (m_unFlags &= ~rfConsuming); };

    /**
     * Set this delay type to non-consuming.
     * 
     * @param consuming new value.
     */
  inline void setDelayNonConsuming(bool nonconsuming) { setDelayConsuming(!nonconsuming); };

    /**
     * Query whether this delay type is consuming.
     * 
     * @return @true if property is set, @false otherwise.
     */
  inline bool isSetDelayConsuming() const { return (m_unFlags & rfConsuming); };

    /**
     * Query whether this delay type is non-consuming.
     * 
     * @return @true if property is set, @false otherwise.
     */
  inline UINTEGER isSetDelayNonConsuming() const { return !isSetDelayConsuming(); };

//     /**
//      * Get a string represantation of this object.
//      * 
//      * @return string representing this object.
//      */
//     STRING toString() const
//     {
//       if(strName.length() > 0)
//         return strName + " [" + strId + "]";
//       else
//         return "[" + strId + "]";
//     }
  };

} } } // close namespaces detail, datamodel & pssalib

#include "Model.h"
#include "SpeciesReference.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /*
   * Get a reactant species reference from the respective list.
   * 
   * Implementation.
   */
 inline const SpeciesReference * Reaction::getReactantsListAt(UINTEGER n) const
  {
    if(n >= m_unReactants)
      return NULL;
    return &m_arSpeciesRefs[n];
  }

  /*
   * Get a reactant species reference from the respective list.
   * 
   * Implementation.
   */
 inline SpeciesReference * Reaction::getReactantsListAt(UINTEGER n)
  {
    if(n >= m_unReactants)
      return NULL;
    return &m_arSpeciesRefs[n];
  }

  /*
   * Get a product species reference from the respective list.
   * 
   * Implementation.
   */
 inline const SpeciesReference * Reaction::getProductsListAt(UINTEGER n) const
  {
    if(n >= (m_unSpeciesRefsCount - m_unReactants))
      return NULL;
    return &m_arSpeciesRefs[m_unReactants + n];
  }

  /*
   * Get a product species reference from the respective list.
   * 
   * Implementation.
   */
 inline SpeciesReference * Reaction::getProductsListAt(UINTEGER n)
  {
    if(n >= (m_unSpeciesRefsCount - m_unReactants))
      return NULL;
    return &m_arSpeciesRefs[m_unReactants + n];
  }

  /*
   * Get a species reference from the list.
   * 
   * Implementation.
   */
 inline const SpeciesReference * Reaction::getSpeciesReferenceAt(UINTEGER n) const
  {
    if(n >= m_unSpeciesRefsCount)
      return NULL;
    return &m_arSpeciesRefs[n];
  }

  /*
   * Get a species reference from the list.
   * 
   * Implementation.
   */
 inline SpeciesReference * Reaction::getSpeciesReferenceAt(UINTEGER n)
  {
    if(n >= m_unSpeciesRefsCount)
      return NULL;
    return &m_arSpeciesRefs[n];
  }

  /*
   * Get the species reference list.
   * 
   * Implementation.
   */
 inline const SpeciesReference * Reaction::getSpeciesReferences() const
  {
    return m_arSpeciesRefs;
  }

  /**
    * Remove the species reference at the given positions in the list.
    * 
    * Implementation.
    */
  inline bool Reaction::removeSpeciesReferencesAt(UINTEGER n)
  {
    if(n >= m_unSpeciesRefsCount)
      return false;
    if(m_unSpeciesRefsCount > (n + 1))
      for(UINTEGER si = n + 1; si < m_unSpeciesRefsCount; ++si)
        m_arSpeciesRefs[si-1] = m_arSpeciesRefs[si];
//       memmove((void*)(m_arSpeciesRefs + n), (void*)(m_arSpeciesRefs + n + 1), (m_unSpeciesRefsCount - n - 1) * sizeof(SpeciesReference));
    --m_unSpeciesRefsCount;
    if(n < m_unReactants)
      --m_unReactants;
    return true;
  }

  /*
   * Swap the species reference at the given positions in the list.
   * 
   * Implementation.
   */
 inline bool Reaction::swapSpeciesReferencesAt(UINTEGER n1, UINTEGER n2)
  {
    if((n1 >= m_unSpeciesRefsCount)||(n2 >= m_unSpeciesRefsCount))
      return false;
    std::swap(m_arSpeciesRefs[n1], m_arSpeciesRefs[n2]);
    return true;
  }

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_REACTION_H_ */
