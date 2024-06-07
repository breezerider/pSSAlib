/**
 * @file DataModel_DM.h
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
 * Gillespie's Direct Method
 */

#ifndef PSSALIB_DATAMODEL_DATAMODEL_DM_H_
#define PSSALIB_DATAMODEL_DATAMODEL_DM_H_

#include "./DataModel.h"
#include "./detail/Subvolume_DM.hpp"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class DataModel_DM
   * @brief Defines the datastructures for the Gillespie's Direct Method.
   *
   * @copydoc DataModel
   */

  class DataModel_DM : public DataModel
  {
  /////////////////////////////////////
  // Constructors
  public:
    // Default constructor
    DataModel_DM();

    //! Copy constructor
    DataModel_DM(DataModel &) = delete;

    // Destructor
  virtual ~DataModel_DM();

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
        new detail::Subvolume_DM);
    };

  /////////////////////////////////////
  // Methods
  public:
    // Subvolumes
    //

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    const detail::Subvolume_DM & getSubvolume(UINTEGER unSubvolumeIdx) const
    {
      return const_cast<const detail::Subvolume_DM &>(
        const_cast<DataModel_DM *>(this)->getSubvolume(unSubvolumeIdx));
    };

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    detail::Subvolume_DM & getSubvolume(UINTEGER unSubvolumeIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unSubvolumeIdx >= m_unSubvolumes)
        throw std::runtime_error("DataModel_DM::getSubvolume() - invalid arguments.");
#endif
      return static_cast<detail::Subvolume_DM &>(*(m_arSubvolumes[unSubvolumeIdx]));
    };

    //! Assignement operator
    DataModel_DM& operator= (const DataModel_DM&) = delete;
  };
}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_DATAMODEL_DM_H_ */
