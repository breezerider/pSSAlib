/**
 * @file DataModel_PDM.h
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
 * Partial Propensity Direct Method (Ramaswamy, 2009)
 */

#ifndef PSSALIB_DATAMODEL_DATAMODEL_PDM_H_
#define PSSALIB_DATAMODEL_DATAMODEL_PDM_H_

#include "./DataModel.h"
#include "./detail/JaggedMatrix.hpp"
#include "./detail/Subvolume_PDM.hpp"
#include "./detail/ReactionWrapper.hpp"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class DataModel_PDM
   * @brief Defines the datastructures for the Partial Propensity Direct Method.
   *
   * @copydoc DataModel
   */
  class DataModel_PDM : public DataModel
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! @internal Struct for storing the row & column of the dependent propensity
    typedef struct tagPropensityIndex
    {
      UINTEGER i, j;

      REAL rate;
      UINTEGER stoichiometry;

      // Constructor
      tagPropensityIndex(UINTEGER _i,
                         UINTEGER _j)
        : i(_i)
        , j(_j)
      {
        // Do nothing
      }

      // Default constructor
      tagPropensityIndex()
        : i(0)
        , j(0)
        , rate(0.0)
        , stoichiometry(0)
      {
        // Do nothing
      }

      // Destructor
      ~tagPropensityIndex()
      {
        // Do nothing
      }

      // Comparison operator<
      inline bool operator<(const tagPropensityIndex& right) const
      {
        if(i == right.i)
          return (bool)(j < right.j);
        return (bool)(i < right.i);
      };

      // operator<< for console output
      friend std::ostream & operator<<(std::ostream &output,
                                       const tagPropensityIndex &pi)
      {
        output << '(' << pi.i << ',' << pi.j << ')';
        return output;
      };
    } PropensityIndex;

  /////////////////////////////////////
  // Constructors
  public:

    //! Default constructor
    DataModel_PDM();

    //! Copy constructor
    DataModel_PDM(DataModel &);

    //! Destructor
  virtual ~DataModel_PDM();

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
        new detail::Subvolume_PDM);
    };

    /**
     * @copydoc DataModel::freeSubvolume()
     */
  virtual void freeSubvolume(detail::Subvolume * sv)
    {
      delete static_cast<detail::Subvolume_PDM *>(sv);
    };

  /////////////////////////////////////
  // Methods
  public:

    /**
     * Clear global data structures.
     */
  virtual void clearStructures()
    {
      arU3.clear();
      aruL.clear();

      // call base class method
      DataModel::clearStructures();
    };

    // Subvolumes
    //

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    const detail::Subvolume_PDM & getSubvolume(UINTEGER unSubvolumeIdx) const
    {
      return const_cast<const detail::Subvolume_PDM &>(
        const_cast<DataModel_PDM *>(this)->getSubvolume(unSubvolumeIdx));
    };

    /**
     * @copydoc DataModel::getSubvolume(UINTEGER)
     */
    detail::Subvolume_PDM & getSubvolume(UINTEGER unSubvolumeIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unSubvolumeIdx >= m_unSubvolumes)
        throw std::runtime_error("DataModel_PDM::getSubvolume() - invalid arguments.");
#endif
      return static_cast<detail::Subvolume_PDM &>(*(m_arSubvolumes[unSubvolumeIdx]));
    };

  ////////////////////////////////
  // Attributes
  public:
    //! Indices of the propensities that need to be updated after
    //! a given reaction has fired.
    detail::JaggedMatrix<PropensityIndex>  arU3;

    //! Look-up table to translate from position in the partial
    //! propensity matrix to reaction index.
    detail::JaggedMatrix<detail::ReactionWrapper *> aruL;
  };

}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_DATAMODEL_PDM_H_ */
