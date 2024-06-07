/**
 * @file DataModel_PSSACR.h
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
 * Declares a datatype containing all the datastructures used by the 
 * Partial Propensity Direct Method with Composition-Rejection sampling
 * (Ramaswamy, 2010)
 */

#ifndef PSSALIB_DATAMODEL_DATAMODEL_PSSACR_H_
#define PSSALIB_DATAMODEL_DATAMODEL_PSSACR_H_

#include "./DataModel_PDM.h"
#include "./CompositionRejectionSamplerData.h"
#include "./detail/Subvolume_PSSACR.hpp"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class DataModel_PSSACR
   * @brief Defines the datastructures for the Partial Propensity Direct Method
   * with Composition-Rejection sampling.
   *
   * @copydoc DataModel
   */
  class DataModel_PSSACR : public DataModel_PDM
  {
  /////////////////////////////////////
  // Constructors
  public:

    // Default constructor
    DataModel_PSSACR();

    //! Copy constructor
    DataModel_PSSACR (const DataModel_PSSACR&) = delete;

    // Destructor
  virtual ~DataModel_PSSACR();

  /////////////////////////////////////
  // Methods
  protected:
    // Subvolumes
    //

    /**
     * @copydoc DataModel::allocateSubvolums
     */
  virtual detail::Subvolume * allocateSubvolume()
    {
      return static_cast<detail::Subvolume *>(
        new detail::Subvolume_PSSACR);
    };

    /**
     * @copydoc DataModel::freeSubvolume()
     */
  virtual void freeSubvolume(detail::Subvolume * sv)
    {
      delete static_cast<detail::Subvolume_PSSACR *>(sv);
    };

  /////////////////////////////////////
  // Methods
  public:
    // Subvolumes
    //

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    const detail::Subvolume_PSSACR & getSubvolume(UINTEGER unSubvolumeIdx) const
    {
      return const_cast<const detail::Subvolume_PSSACR &>(
        const_cast<DataModel_PSSACR *>(this)->getSubvolume(unSubvolumeIdx));
    };

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    detail::Subvolume_PSSACR & getSubvolume(UINTEGER unSubvolumeIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unSubvolumeIdx >= m_unSubvolumes)
        throw std::runtime_error("DataModel_PSSACR::getSubvolume() - invalid arguments.");
#endif
      return static_cast<detail::Subvolume_PSSACR &>(*(m_arSubvolumes[unSubvolumeIdx]));
    };

    //! Assignement operator
    DataModel_PSSACR& operator= (const DataModel_PSSACR&) = delete;

  };
}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_DATAMODEL_PSSACR_H_ */
