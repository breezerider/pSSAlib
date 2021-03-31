/**
 * @file CompositionRejectionSamplerData.h
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
 * Sampler data for the Composition-Rejection method
 */

#ifndef PSSALIB_DATAMODEL_COMPOSITION_REJECTION_SAMPLER_DATA_H_
#define PSSALIB_DATAMODEL_COMPOSITION_REJECTION_SAMPLER_DATA_H_

#include "../typedefs.h"
#include "./PSSACR_Bins.h"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class CompositionRejectionSamplerData
   * @brief Stores binned dat for composition-rejection sampling.
   */
  class CompositionRejectionSamplerData
  {
  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    CompositionRejectionSamplerData()
      : minValue(0.0)
    {
      // Do nothing
    }

    //! Copy constructor
    CompositionRejectionSamplerData(
      const CompositionRejectionSamplerData &other)
      : minValue(other.minValue)
      , bins(other.bins)
    {
      // Do nothing
    }

  /////////////////////////////////////
  // Methods
  public:
    //! assignment operator
  inline CompositionRejectionSamplerData &operator=(
      const CompositionRejectionSamplerData &other)
    {
      minValue = other.minValue;
      bins = const_cast<PSSACR_Bins&>(other.bins);
      return *this;
    }

    /**
     * Update value
     * 
     * @param bin_no_new
     * @param idx
     * @param val
     */
  inline void updateValue(UINTEGER bin_no_new, 
                            UINTEGER idx, REAL val)
    {
      bins.updateValue(bin_no_new, idx, val);
    }

    /**
     * Clear binned values
     */
  inline void clear()
    {
      minValue = 0.0;
      bins.clear();
    }

  ////////////////////////////////
  // Attributes
  public:
    REAL        minValue;
    PSSACR_Bins bins;
  };

}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_COMPOSITION_REJECTION_SAMPLER_DATA_H_ */
