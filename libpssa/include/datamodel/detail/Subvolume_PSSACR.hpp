/**
 * @file Subvolume_PSSACR.hpp
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
 * Declares a container for subvolume variables required by the
 * PSSA with Composition-Rejection Sampling (Ramaswamy, 2010)
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PSSACR_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PSSACR_HPP_

#include "../../typedefs.h"
#include "Subvolume_PDM.hpp"

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel_PSSACR;
  
namespace detail
{
  /**
   * @copydoc Subvolume
   */
  class Subvolume_PSSACR : public Subvolume_PDM
  {
  ////////////////////////////////
  // Friends
  public:
    friend class DataModel_PSSACR;

  ////////////////////////////////
  // Attributes
  protected:
    // Reactions
    //

    // 
    CompositionRejectionSamplerData *m_crsdPi;

  ////////////////////////////////
  // Attributes
  public:
    // Reactions
    //

    //
    CompositionRejectionSamplerData crsdSigma;

  ////////////////////////////////
  // Constructors
  public:
    //! Constructor
    Subvolume_PSSACR()
      : m_crsdPi(NULL)
    {
      // Do nothing
    }

    //! Destructor
  virtual ~Subvolume_PSSACR()
    {
      // Clean-up
      free_PSSACR();
    }

  ////////////////////////////////
  // Methods
  private:

    /*
     * Free memory
     */
    void free_PSSACR()
    {
      if(NULL != m_crsdPi)
      {
        delete [] m_crsdPi;
        m_crsdPi = NULL;
      }
    };

  ////////////////////////////////
  // Methods
  protected:
    /**
     * Reset all properties' values.
     */
  virtual void free()
    {
      // free memory
      free_PSSACR();

      // call base class method
      Subvolume_PDM::free();
    };

    /**
     * @copydoc Subvolume::allocate(UINTEGER,UINTEGER,BYTE)
     */
  virtual void allocate(UINTEGER reactions, UINTEGER species, BYTE dims)
    {
      // call base class method
      Subvolume_PDM::allocate(reactions, species, dims);

      // allocate memory
      m_crsdPi = new CompositionRejectionSamplerData[species];
    };

    /**
     * @copydoc Subvolume::clear(UINTEGER,UINTEGER)
     */
  virtual void clear(UINTEGER reactions, UINTEGER species)
    {
      for(UINTEGER si = 0; si < species; ++si)
        m_crsdPi[si].clear();
      crsdSigma.clear();

      // call base class method
      Subvolume_PDM::clear(reactions, species);
    };

  ////////////////////////////////
  // Methods
  public:

    /**
     * Get partial propensity bins for species 
     * 
     * @param index Species index in the model
     * @return Partial propensity bins
     */
  inline CompositionRejectionSamplerData & crsdPi(UINTEGER index)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= unSpecies)
        throw std::runtime_error("Subvolume_PSSACR::crsdPi() - invalid arguments.");
#endif
      return m_crsdPi[index];
    }
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PSSACR_HPP_ */
