/**
 * @file DataModel_SPDM.h
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
 * Declares a container for all the data structures required by the 
 * Sorting Partial Propensity Direct Method (Ramaswamy, 2009)
 */

#ifndef PSSALIB_DATAMODEL_DATAMODEL_SPDM_H_
#define PSSALIB_DATAMODEL_DATAMODEL_SPDM_H_

#include "./DataModel_PDM.h"
#include "./detail/Subvolume_SPDM.hpp"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class DataModel_SPDM
   * @brief Defines the datastructures for the Sorting Partial Propensity 
   * Direct Method.
   *
   * @copydoc DataModel
   */
  class DataModel_SPDM : public DataModel_PDM
  {
  /////////////////////////////////////
  // Constructors
  public:

    // Default constructor
    DataModel_SPDM();

    // Copy constructor
    DataModel_SPDM(DataModel &);

    // Destructor
  virtual ~DataModel_SPDM();

  /////////////////////////////////////
  // Methods
  protected:
    // Subvolumes
    //

    /**
     * @copydoc DataModel::allocateSubvolume()
     */
  virtual detail::Subvolume * allocateSubvolume()
    {
      return static_cast<detail::Subvolume *>(
        new detail::Subvolume_SPDM);//internalAllocateSubvolumes<detail::Subvolume_SPDM>(n));
    };

    /**
     * @copydoc DataModel::freeSubvolume()
     */
  virtual void freeSubvolume(detail::Subvolume * sv)
    {
      delete static_cast<detail::Subvolume_SPDM *>(sv);
    };

  /////////////////////////////////////
  // Methods
  public:
    // Subvolumes
    //

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    const detail::Subvolume_SPDM & getSubvolume(UINTEGER unSubvolumeIdx) const
    {
      return const_cast<const detail::Subvolume_SPDM &>(
        const_cast<DataModel_SPDM *>(this)->getSubvolume(unSubvolumeIdx));
    };

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    detail::Subvolume_SPDM & getSubvolume(UINTEGER unSubvolumeIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unSubvolumeIdx >= m_unSubvolumes)
        throw std::runtime_error("DataModel_SPDM::getSubvolume() - invalid arguments.");
#endif
      return static_cast<detail::Subvolume_SPDM &>(*(m_arSubvolumes[unSubvolumeIdx]));
    };

  ////////////////////////////////
  // Attributes
  public:
    // Reactions
    //

    // Position
    std::size_t rowIndex, colIndex;
  };
}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_DATAMODEL_SPDM_H_ */
