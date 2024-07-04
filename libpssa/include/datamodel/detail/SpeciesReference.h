/**
 * @file SpeciesReference.h
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
 * Declares a container for properties of SBML SpeciesReference objects that are
 * relevant for simulations using Stochastic Simulation Algorithms
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SPECIESREFERENCE_H_
#define PSSALIB_DATAMODEL_DETAIL_SPECIESREFERENCE_H_

#include "Base.hpp"

#ifndef PSSALIB_SPECIES_ID_RESERVOIR
#define PSSALIB_SPECIES_ID_RESERVOIR std::numeric_limits<UINTEGER>::max()
#endif

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  // Forward declarations
  class Reaction;

  /**
   * @class SpeciesReference
   * @brief Declares a datatype for a refernce to chemical species
   */
  class SpeciesReference : public Base
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagSpeciesReferenceFlags
    {
      srfConstant = bfAll + 1,

      srfAll = srfConstant
    } SpeciesReferenceFlags;

  ////////////////////////////////
  // Attributes
  protected:

    //! species index
    UINTEGER m_unSpeciesIndex;

    //! species stoichiometry
    BYTE m_ucSpeciesStoichiometry;

    //! pointer to parent reaction
    Reaction * m_ptrReaction;

  /////////////////////////////////////
  // Constructors
  public:
    // Default constructor
    SpeciesReference(Reaction * reaction = NULL)
      : Base()
      , m_unSpeciesIndex(0)
      , m_ucSpeciesStoichiometry(0)
      , m_ptrReaction(reaction)
    {
      // Do nothing
    };
#ifdef HAVE_LIBSBML
    // Default constructor
    SpeciesReference(const LIBSBML_CPP_NAMESPACE::SpeciesReference * speciesReference, Reaction * reaction)
      : Base()
      , m_unSpeciesIndex(0)
      , m_ucSpeciesStoichiometry(0)
      , m_ptrReaction(NULL)
    {
      assign(speciesReference, reaction);
    };
#endif
    // Copy constructor
    SpeciesReference(const SpeciesReference & right)
      : Base(right)
      , m_unSpeciesIndex(right.m_unSpeciesIndex)
      , m_ucSpeciesStoichiometry(right.m_ucSpeciesStoichiometry)
      , m_ptrReaction(NULL)
    {
      // Do nothing
    };

    // Destructor
    ~SpeciesReference()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  public:
    /**
     * Reset all properties' values.
     */
  virtual void unset()
    {
      // call the base class method
      Base::unset();

      m_unSpeciesIndex = 0;
      m_ucSpeciesStoichiometry = 0;
    };
#ifdef HAVE_LIBSBML
    /**
     * Assign properties from an SBML SpeciesReference object.
     * 
     * @param speciesReference an SBML SpeciesReference object.
     * @param reaction @c Reaction container object.
     * @return @true if assignment was successful, @false otherwise.
     */
    bool assign(const LIBSBML_CPP_NAMESPACE::SpeciesReference * speciesReference, Reaction * reaction);
#endif
    /**
     * Get the species reference stoichiometry according to its role.
     * 
     * @return current value, always positive.
     */
  inline INTEGER getStoichiometry() const
    {
        return ((INTEGER)m_ucSpeciesStoichiometry);
    };

    /**
     * Get the species reference stoichiometry as an absolute value.
     * 
     * @return current value.
     */
  inline UINTEGER getStoichiometryAbs() const
    {
      return ((UINTEGER)m_ucSpeciesStoichiometry);
    };

    /**
     * Set the species stoichiometry.
     * 
     * @param stoichiometry new value.
     */
  inline void setStoichiometry(const BYTE stoichiometry)
    {
      m_ucSpeciesStoichiometry = stoichiometry;
    };

    /**
     * Get the species index.
     * 
     * @return current value.
     */
  inline UINTEGER getIndex() const
    {
      return m_unSpeciesIndex;
    };

    /**
     * Set the species index.
     * 
     * @param index new value.
     */
  inline void setIndex(const UINTEGER index)
    {
      m_unSpeciesIndex = index;
    };

    /**
     * Get the pointer to parent reaction.
     *
     * @return current value.
     */
  inline const Reaction * getReaction() const
    {
      return m_ptrReaction;
    };

    /**
     * Set the pointer to parent reaction.
     *
     * @param reaction new value.
     */
  inline void setReaction(Reaction * const reaction)
    {
      m_ptrReaction = reaction;
    };

    /**
     * Query whether this species reference points to a constant species.
     * 
     * @return @true if target species is set to be kept constant, @false otherwise.
     */
  inline bool isConstant() const { return (m_unFlags & srfConstant); };

    /**
     * Set whether this species reference points to a constant species.
     * 
     * @param c @true if target species is to be kept constant, @false otherwise.
     */
  inline void setConstant(const bool c)
    {
      c ? (m_unFlags |= srfConstant) : (m_unFlags &= ~srfConstant);
    };

    /**
     * @internal Create a reference to reservoid species
     */
  inline void makeReservoir()
    {
      m_unFlags |= srfConstant; m_unSpeciesIndex = PSSALIB_SPECIES_ID_RESERVOIR;
    }

    /**
     * @internal Create a reference to reservoid species
     * 
     * @return @true if target species is reservoir species, @false otherwise.
     */
  inline bool isReservoir() const
    {
      return (m_unFlags & srfConstant) && (PSSALIB_SPECIES_ID_RESERVOIR == m_unSpeciesIndex);
    }

    /**
     * @internal Used to match this instance against another one.
     * 
     * @param o Anoother instance of this class.
     */
    bool operator()(const SpeciesReference & o) const
    {
      return (o.m_unSpeciesIndex == m_unSpeciesIndex);
    };

    /**
     * Get a string represantation of this object.
     * 
     * @return string representing this object.
     */
    STRING toString() const;

  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SPECIESREFERENCE_H_ */
