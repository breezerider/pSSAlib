/**
 * @file ReactionWrapper.hpp
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
 * Wrapper over a Rection or Species object that provides a unified interface
 * for treating ordinary, reversible and diffusion reactions
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_REACTIONWRAPPER_HPP_
#define PSSALIB_DATAMODEL_DETAIL_REACTIONWRAPPER_HPP_

#include "Reaction.h"
#include "Species.h"
#include "SpeciesReference.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /**
   * @class ReactionWrapper
   * @brief Decorator for the Reaction datatype
   */
  class ReactionWrapper
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    typedef enum tagReactionWrapperFlags
    {
      rwfReverse = 0x01,
      rwfDelayedNonCosumingUpdate = 0x02,
      rwfDiffusion = 0x4,

      rwfAll = rwfReverse + rwfDelayedNonCosumingUpdate + rwfDiffusion
    } ReactionWrapperFlags;

    //! Pointer to underlying type
    typedef union tagReactionWrapperComponent
    {
      Reaction * ptrReaction;
      Species  * ptrSpecies;
    } ReactionWrapperComponent;

  ////////////////////////////////
  // Attributes
  protected:
    //! component
    ReactionWrapperComponent component;

    //! serial number
    UINTEGER unSerialNumber;

    //! flags
    UINTEGER unFlags;

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
  explicit ReactionWrapper(UINTEGER serialNumber)
      : unSerialNumber(serialNumber)
      , unFlags(0)
    {
      // Do nothing
    };

    //! Constructor (ordinary reaction)
    ReactionWrapper(Reaction * reaction, UINTEGER serialNumber, bool reverse = false)
      : unSerialNumber(serialNumber)
      , unFlags(reverse ? rwfReverse : 0)
    {
      // Do nothing
      component.ptrReaction = reaction;
    };

    //! Constructor (diffusion reaction)
    ReactionWrapper(Species * species, UINTEGER serialNumber)
      : unSerialNumber(serialNumber)
      , unFlags(rwfDiffusion)
    {
      // Do nothing
      component.ptrSpecies = species;
    };

    //! Copy constructor
    ReactionWrapper(const ReactionWrapper & right)
      : unSerialNumber(right.unSerialNumber)
      , unFlags(right.unFlags)
    {
      if(unFlags & rwfDiffusion)
        component.ptrSpecies = right.component.ptrSpecies;
      else
        component.ptrReaction = right.component.ptrReaction;
    };

    //! Destructor
    ~ReactionWrapper()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  public:
    /**
     * Query whether this reaction decorator is for a reverse reaction.
     * 
     * @return @true if this is a reverse reaction decorator, @false otherwise.
     */
  inline bool isReverse() const { return (unFlags & rwfReverse); };

    /**
     * Query whether this is a reaction decorator for diffusion reaction.
     * 
     * @return @true if this is a diffusion reaction decorator, @false otherwise.
     */
  inline bool isDiffusive() const { return (unFlags & rwfDiffusion); };

    /**
     * Query whether this reaction requires a delayed update.
     * 
     * @return @true if a delayed update is required, @false otherwise.
     */
  inline bool isSetDelay() const { return ((unFlags & rwfDiffusion) ? false : component.ptrReaction->isSetDelay()); };

    /**
     * Query whether this reaction's delay is consuming.
     * 
     * @return @true if the delay is consuming, @false otherwise.
     */
  inline bool isSetDelayConsuming() const { return ((unFlags & rwfDiffusion) ? false : component.ptrReaction->isSetDelayConsuming()); };

    /**
     * Get the reaction delay.
     * 
     * @return current value.
     */
  inline REAL getDelay() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support delays");
      else
#endif
        return component.ptrReaction->getDelay();
    };

    /**
     * Get the reaction rate.
     * 
     * @return current value.
     */
  inline REAL getRate() const
    {
      if(unFlags & rwfDiffusion)
        return component.ptrSpecies->getDiffusionConstant();
      else if(unFlags & rwfReverse)
        return component.ptrReaction->getReverseRate();
      else
        return component.ptrReaction->getForwardRate();
    };

    /**
     * Get the reaction serial number.
     * 
     * @return current number
     */
  inline UINTEGER getSerialNumber() const
    {
      return unSerialNumber;
    };

   /**
     * @copydoc Reaction::getReactantsCount()
     */
  inline UINTEGER getReactantsCount() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getReactantsCount()");
      else
#endif
           if(unFlags & rwfReverse)
        return component.ptrReaction->getProductsCount();
      else
        return component.ptrReaction->getReactantsCount();
    };

    /**
     * @copydoc Reaction::getProductsCount()
     */
  inline UINTEGER getProductsCount() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getProductsCount()");
      else
#endif
           if(unFlags & rwfReverse)
        return component.ptrReaction->getReactantsCount();
      else
        return component.ptrReaction->getProductsCount();
    };

    /**
     * @copydoc Reaction::getSpeciesReferencesCount()
     */
  inline UINTEGER getSpeciesReferencesCount() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getSpeciesReferencesCount()");
      else
#endif
        return component.ptrReaction->getSpeciesReferencesCount();
    };

    /**
     * @copydoc Reaction::getReactantsListAt(UINTEGER)
     */
  inline const SpeciesReference * getReactantsListAt(UINTEGER n) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getReactantsListAt()");
      else
#endif
           if(unFlags & rwfReverse)
        return component.ptrReaction->getProductsListAt(n);
      else
        return component.ptrReaction->getReactantsListAt(n);
    };

    /**
     * @copydoc Reaction::getProductsListAt(UINTEGER)
     */
  inline const SpeciesReference * getProductsListAt(UINTEGER n) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getProductsListAt()");
      else
#endif
           if(unFlags & rwfReverse)
        return component.ptrReaction->getReactantsListAt(n);
      else
        return component.ptrReaction->getProductsListAt(n);
    };

    /**
     * @copydoc Reaction::getSpeciesReferenceAt(UINTEGER)
     */
  inline const SpeciesReference * getSpeciesReferenceAt(UINTEGER n) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getSpeciesReferenceAt()");
      else
#endif
           if(unFlags & rwfReverse)
        return component.ptrReaction->getSpeciesReferenceAt(
          component.ptrReaction->getSpeciesReferencesCount() - ++n);
      else
        return component.ptrReaction->getSpeciesReferenceAt(n);
    };

    /**
     * @copydoc Reaction::getSpeciesReferences()
     */
  inline const SpeciesReference * getSpeciesReferences() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support getSpeciesReferences()");
      else
#endif
        return component.ptrReaction->getSpeciesReferences();
    }

    /**
     * @copydoc Reaction::swapSpeciesReferencesAt(UINTEGER,UINTEGER)
     */
  inline bool swapSpeciesReferencesAt(UINTEGER n1, UINTEGER n2)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unFlags & rwfDiffusion)
        throw std::runtime_error("diffusion reactions do not support swapSpeciesReferencesAt()");
      else
#endif
        return component.ptrReaction->swapSpeciesReferencesAt(n1, n2);
    }

    /**
     * Get the species for a diffusion reaction.
     * 
     * @param n Index of the species reference to retrieve.
     * @return SpeciesReference object at given position in the
     * list or NULL if the index exceeds the upper bound.
     */
  inline const Species * getSpecies() const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(0 == (unFlags & rwfDiffusion))
        throw std::runtime_error("reactions do not support getSpecies()");
      else
#endif
        return component.ptrSpecies;
    };

    void getSymbolicRepresentation(std::ostream & os) const
    {
      if(unFlags & rwfDiffusion)
      {
        os << component.ptrSpecies->toString() << " --" << component.ptrSpecies->getDiffusionConstant() << "--> " << component.ptrSpecies->toString();
      }
      else
      {
        for (UINTEGER ri = 0 ; ri < getReactantsCount(); ri++)
        {
          const SpeciesReference *spRef = getReactantsListAt(ri);
          if(spRef->isReservoir())
            os << "[]";
          else
          {
//             if(spRef->getStoichiometryAbs() > 1)
              os << spRef->getStoichiometryAbs() << " * ";
            os << component.ptrReaction->getModel()->getSpecies(spRef->getIndex())->toString();
          }
          os << ((ri < (getReactantsCount() - 1)) ? " + " : " ");
        }
        os << "--" << getRate() << "--> ";
        for (UINTEGER pi = 0 ; pi < getProductsCount(); pi++)
        {
          const SpeciesReference *spRef = getProductsListAt(pi);
          if(spRef->isReservoir())
            os << "[]";
          else
          {
//             if(spRef->getStoichiometryAbs() > 1)
              os << spRef->getStoichiometryAbs() << " * ";
            os << component.ptrReaction->getModel()->getSpecies(spRef->getIndex())->toString();
          }
          os << ((pi < (getProductsCount() - 1)) ? " + " : " ");
        }
      }
    };

    /**
     * Get a string represantation of this object.
     * 
     * @return string representing this object.
     */
    STRING toString() const
    {
      if(unFlags & rwfDiffusion)
      {
        return component.ptrSpecies->toString();
      }
      else
      {
        STRINGSTREAM ssTemp;
        getSymbolicRepresentation(ssTemp);
        return ssTemp.str();
      }
    };
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_REACTIONWRAPPER_HPP_ */
