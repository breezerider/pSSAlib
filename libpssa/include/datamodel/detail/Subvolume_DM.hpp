/**
 * @file Subvolume_DM.hpp
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
 * Declares a container for subvolume variables used by the
 * Gillespie's Direct Method
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_DM_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_DM_HPP_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel_DM;

namespace detail
{
  /**
   * @class Subvolume
   * @brief Container defining a subreactor state
   */
  class Subvolume_DM : public Subvolume
  {
  ////////////////////////////////
  // Friends
  public:
    friend class DataModel_DM;

  ////////////////////////////////
  // Attributes
  protected:
    // Reactions
    //

    //! Propensities
    REAL *ardPi;

  ////////////////////////////////
  // Constructors
  public:
    //! Constructor
    Subvolume_DM()
      : ardPi(NULL)
    {
      // Do nothing
    }

    //! Destructor
  virtual ~Subvolume_DM()
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
      if(NULL != ardPi)
      {
        delete [] ardPi;
        ardPi = NULL;
      }

      Subvolume::free();
    };

    /**
     * Allocate the subvolume datastructure
     * 
     * @param reactions number of reactions in the model
     * @param species number of species in the model
     * @param dims number of spatial dimensions
     */
  virtual void allocate(UINTEGER reactions, UINTEGER species, BYTE dims)
    {
      // call base class method
      Subvolume::allocate(reactions, species, dims);

      // allocate memory
      ardPi = new REAL[reactions];
      std::fill_n(ardPi, reactions, REAL(0.0));
    };

  ////////////////////////////////
  // Methods
  public:

    /**
     * Get reaction propensity
     * 
     * @param index Reaction index in the model
     * @return Reaction propensity
     */
  inline REAL & propensity(UINTEGER index)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= unReactions)
        throw std::runtime_error("DataModel::getSubvolume() - invalid arguments.");
#endif
      return ardPi[index];
    }

  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_DM_HPP_ */
