/**
 * @file Subvolume.hpp
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
 * Declares a container for common subvolume variables
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_HPP_

#include "../../typedefs.h"
#include "SpeciesReference.h"

#ifndef PSSALIB_NO_BOUNDS_CHECKS
#if __GNUC__ > 4 || \
    (__GNUC__ == 4 && (__GNUC_MINOR__ >= 2))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#endif

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel;

namespace detail
{
  /**
   * @class Subvolume
   * @brief Container defining a subreactor state
   */
  class Subvolume
  {
  ////////////////////////////////
  // Friends
  public:
    friend class pssalib::datamodel::DataModel;

  ////////////////////////////////
  // Attributes
  protected:
    // Species
    //

    //! Vector of current species population
    UINTEGER *arunPopulation;

    // Subvolume
    //

    //! Indexes of the neighboring subvolumes
    UINTEGER *arNeighbouringSubvolumes;

    // Debugging
#ifndef PSSALIB_NO_BOUNDS_CHECKS
    UINTEGER unSpecies, unReactions;
    BYTE     uDims;
#endif

  ////////////////////////////////
  // Attributes
  public:
    // Reactions
    //

    //! Total propensity
    REAL     dTotalPropensity;

  ////////////////////////////////
  // Constructors
  public:
    //! Constructor
    Subvolume()
      : arunPopulation(NULL)
      , arNeighbouringSubvolumes(NULL)
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      , unSpecies(0)
      , unReactions(0)
      , uDims(0)
#endif
      , dTotalPropensity(0.0)
    {
      // Do nothing
    }

    //! Destructor
  virtual ~Subvolume()
    {
      // Clean-up
      free();
    }

  ////////////////////////////////
  // Methods
  protected:
    /**
     * Reset all properties' values.
     */
  virtual void free()
    {
      if(NULL != arunPopulation)
      {
        delete [] arunPopulation;
        arunPopulation = NULL;
      }
      if(NULL != arNeighbouringSubvolumes)
      {
        delete [] arNeighbouringSubvolumes;
        arNeighbouringSubvolumes = NULL;
      }
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      unSpecies = 0;
      unReactions = 0;
      uDims = 0;
#endif
      dTotalPropensity = 0.0;
    };

    /**
     * Allocate memory for subvolume data structures.
     * 
     * @param reactions number of reactions in the model.
     * @param species number of species in the model.
     * @param dims number of spatial dimensions.
     */
  virtual void allocate(UINTEGER reactions, UINTEGER species, BYTE dims)
    {
      // clean up
      free();

      // allocate memory
      arunPopulation = new UINTEGER[species];
      if(0 != dims)
        arNeighbouringSubvolumes = new UINTEGER[2*dims];

#ifndef PSSALIB_NO_BOUNDS_CHECKS
      unReactions = reactions;
      unSpecies = species;
      uDims = dims;
#endif
    };

    /**
     * Clear all simulation states variables.
     * 
     * @param reactions number of reactions in the model.
     * @param species number of species in the model.
     */
  virtual void clear(UINTEGER reactions, UINTEGER species)
    {
      memset(arunPopulation, (unsigned char)0, sizeof(UINTEGER)*species);
    };

  ////////////////////////////////
  // Methods
  public:

    /**
     * Update species population
     * 
     * @param index Species index in the model
     * @param change Change in species population
     */
  inline void population_update(const UINTEGER index, const INTEGER change)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((index >= unSpecies)||
         ((change < 0)&&(arunPopulation[index] < std::abs(change))))
        throw std::runtime_error("Subvolume::population_update() - invalid arguments.");
#endif
      arunPopulation[index] += change;
    }

    /**
     * Update species population
     * 
     * @param sr Pointer to species reference in the model
     * @param grow @true if growing population, @false otherwise
     */
  inline void population_update(const SpeciesReference * sr, const bool grow)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((NULL==sr)||(sr->getIndex() >= unSpecies)||
         ((!grow)&&(arunPopulation[sr->getIndex()] < sr->getStoichiometryAbs())))
        throw std::runtime_error("Subvolume::population_update() - invalid arguments.");
#endif
      if(grow)
        arunPopulation[sr->getIndex()] += sr->getStoichiometryAbs();
      else
        arunPopulation[sr->getIndex()] -= sr->getStoichiometryAbs();
    }

    /**
     * Get species population
     * 
     * @param index Species index in the model
     * @return Species population
     */
  inline UINTEGER population(UINTEGER index) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= unSpecies)
        throw std::runtime_error("Subvolume::population() - invalid arguments.");
#endif
      return arunPopulation[index];
    }

    /**
     * Get subvolume neightbour
     * 
     * @param index Neightbour index
     * @return Neightbouting subvolume index
     */
inline UINTEGER neighbour(UINTEGER index) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= 2*uDims)
        throw std::runtime_error("Subvolume::neighbour() - invalid arguments.");
#endif
      return arNeighbouringSubvolumes[index];
    }

    /**
     * Get a string represantation of this object.
     * 
     * @return string representing this object.
     */
    STRING toString() const
    {
      return STRING("Subvolume");
    };
  };

} } } // close namespaces detail, datamodel & pssalib

#ifndef PSSALIB_NO_BOUNDS_CHECKS
#if __GNUC__ > 4 || \
    (__GNUC__ == 4 && (__GNUC_MINOR__ >= 2))
#pragma GCC diagnostic pop
#endif
#endif

#endif /* PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_HPP_ */

