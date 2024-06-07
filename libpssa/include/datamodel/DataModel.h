/**
 * @file DataModel.h
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
 * Declares a container for common data structures shared by  
 * all Stochastic Simulation Algorithms
 */

#ifndef PSSALIB_DATAMODEL_DATAMODEL_H_
#define PSSALIB_DATAMODEL_DATAMODEL_H_

#include "../typedefs.h"
#include "./CompositionRejectionSamplerData.h"

#include "./detail/Model.h"
#include "./detail/Subvolume.hpp"
#include "./detail/ReactionWrapper.hpp"
#include "./detail/VolumeDecomposition.hpp"

namespace pssalib
{
namespace datamodel
{
  /**
   * @class DataModel
   * @brief This class serves as a base class for the respective 
   * datastructures of the methods offered by this library.
   * 
   * @details Contains a method-specific model definition along with the
   * domain subdivision into homogeneous chemical reactors (subvolumes).
   */
  class DataModel : public detail::Model
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagDataModelFlags
    {
      dmfBCReflexive = mfAll + 1,

      dmfAll = dmfBCReflexive
    } DataModelFlags;

    //! @internal Struct for delayed reactions
    typedef struct tagDelayedReaction
    {
      //! Reaction index
      UINTEGER index;
      //! Time-point when it's scheduled to fire
      REAL     time;

      //! Default constructor
      tagDelayedReaction()
        : index(0)
        , time(0.0)
      {
        // Do nothing
      }

      //! Constructor
      tagDelayedReaction(UINTEGER i, REAL t)
        : index(i)
        , time(t)
      {
        // Do nothing
      }

      //! Destructor
      ~tagDelayedReaction()
      {
        // Do nothing
      }

      //! Comparison operator less than
      inline bool operator<(const tagDelayedReaction &right) const
      {
        return (time < right.time);
      };
    } DelayedReaction;

    //! @internal Struct for handling subreactor boundaries
    struct tagBCHelper
    {
      virtual UINTEGER prev(UINTEGER curr, UINTEGER len) const = 0;
      virtual UINTEGER next(UINTEGER curr, UINTEGER len) const = 0;
    };

    //! @internal Struct for periodic BCs
    struct tagPeriodicBCHelper : public tagBCHelper
    {
      virtual UINTEGER prev(UINTEGER curr, UINTEGER len) const
      {
        return (curr - 1 + len) % len;
      };

      virtual UINTEGER next(UINTEGER curr, UINTEGER len) const
      {
        return (curr + 1 + len) % len;
      };
    };

    //! @internal Struct for reflexive BCs
    struct tagReflexiveBCHelper : public tagBCHelper
    {
      virtual UINTEGER prev(UINTEGER curr, UINTEGER len) const
      {
        return (curr > 0) ? (curr - 1) : (UINTEGER)0;
      };

      virtual UINTEGER next(UINTEGER curr, UINTEGER len) const
      {
        return ((curr + 1) < len) ? (curr + 1) : (len - 1);
      };
    };

  /////////////////////////////////////
  // Constructors
  public:
    // Default constructor
    DataModel();

    //! Copy constructor
    DataModel(DataModel &) = delete;

    // Destructor
  virtual ~DataModel();

  /////////////////////////////////////
  // Methods
  protected:
    // Subvolumes
    //

    /**
     * Allocate a subvolume.
     * 
     * @return Pointer to the allocated object.
     */
  virtual detail::Subvolume * allocateSubvolume()
    {
      return new detail::Subvolume;
      //internalAllocateSubvolumes(n);
    };

    /**
     * Free memory for a subvolume.
     * 
     * @param sv Pointer to the allocated object.
     */
  virtual void freeSubvolume(detail::Subvolume * sv)
    {
      delete sv;
    };

    /**
     * Free the subvolume containers.
     */
    void freeSubvolumes()
    {
      if(NULL != m_arSubvolumes)
      {
        for(UINTEGER svi = 0; svi < m_unSubvolumes; ++svi)
          if(NULL != m_arSubvolumes[svi])
            freeSubvolume(m_arSubvolumes[svi]);
        delete [] m_arSubvolumes;
      }
      m_unSubvolumes = 0;
      m_arSubvolumes = NULL;
    };

  /////////////////////////////////////
  // Methods
  public:
    /**
     * Release the allocated data structures.
     */
  virtual void free();

    /**
     * Clear both global and subvolume data structures.
     */
  virtual void clear()
    {
      clearSubvolumes();
      clearStructures();
    }

    /**
     * Clear subvolume data structures.
     */
  inline void clearSubvolumes()
    {
      for(UINTEGER svi = 0; svi < m_unSubvolumes; ++svi)
      {
        getSubvolume(svi).clear(m_unReactions, m_unSpecies);
        if(0 != m_unSpecies)
          getSubvolume(svi).arunPopulation[0] = 1;
      }
    };

    /**
     * Clear global data structures.
     */
  virtual void clearStructures()
    {
      // Do nothing
    };

    /**
     * Swap member values with another instance
     */
  virtual void swap(DataModel & other)
    {
      Model::swap(other);

      std::swap(m_arReactionWrappers, other.m_arReactionWrappers);
      std::swap(m_unReactionWrappers, other.m_unReactionWrappers);

      std::swap(m_arSubvolumes, other.m_arSubvolumes);
      std::swap(m_unSubvolumes, other.m_unSubvolumes);

      std::swap(m_arunDims, other.m_arunDims);
      std::swap(m_uDims, other.m_uDims);

//       if(m_uDims > 0)
//       setupVolumeDecomposition(m_uDims, m_arunDims, (m_unFlags & dmfBCReflexive) ?
//                                 detail::BC_Reflexive : detail::BC_Periodic);
//       else
//         setupReactorVolume();
    }

//     void setupReactorVolume();

    void setup(BYTE dims, const UINTEGER * pDims, const detail::BoundaryConditionsType & bc);

    void setupVolumeDecomposition(UINTEGER subvolumes, const detail::BoundaryConditionsType & bc);

    void setupReactionWrappers(UINTEGER subvolumes);

  inline void setupPopulation(UINTEGER ** initAmounts)
    {
      for(UINTEGER si = 0; si < m_unSpecies; ++si)
          for(UINTEGER svi = 0; svi < m_unSubvolumes; ++svi)
            getSubvolume(svi).arunPopulation[si] = initAmounts[svi][si];
    }

    // Reactions
    //

    /**
     * Get a reaction wrapper at a given position in the list.
     * 
     * @param unRWIdx List index of the reaction wrapper.
     * @return Reference to the ReactionWrapper object.
     */
  inline const detail::ReactionWrapper & getReactionWrapper(UINTEGER unRWIdx) const
    {
      return const_cast<const detail::ReactionWrapper &>(const_cast<DataModel *>(this)->getReactionWrapper(unRWIdx));
    };

    /**
     * Get a reaction wrapper at a given position in the list.
     * 
     * @param unRWIdx List index of the reaction wrapper.
     * @return Reference to the ReactionWrapper object.
     */
  inline detail::ReactionWrapper & getReactionWrapper(UINTEGER unRWIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unRWIdx >= m_unReactionWrappers)
        throw std::runtime_error("DataModel::getReactionWrapper() - invalid arguments.");
#endif
      return m_arReactionWrappers[unRWIdx];
    };

    /**
     * Get number of reaction wrappers.
     * 
     * @return Total number of ReactionWrapper objects
     * associated with this DataModel.
     */
  inline UINTEGER getReactionWrappersCount() const
    {
      return m_unReactionWrappers;
    };

    // Subvolumes
    //

    /**
     * Get number of spatial dimensions.
     * 
     * @return Number of spatial dimensions used by
     * this DataModel in reactor volume subdivision.
     */
  inline BYTE getDimsCount() const
    {
      return m_uDims;
    };

    /**
     * Get number of subvolumes along a given spatial dimension.
     * 
     * @param d Spatial dimension index (must be less than @see getDimsCount())
     * @return Number of spatial subvolumes along this dimension.
     */
  inline UINTEGER getDims(BYTE d) const
    {
      if(d < m_uDims)
        return m_arunDims[d];
      else
        return 0;
    };

    /**
     * Get all spatial dimensions.
     * 
     * @return Vector of spatial dimensions.
     */
  inline const UINTEGER * getDims() const
    {
      return m_arunDims;
    };

    /**
     * Get a subvolume at a given position in the list.
     * 
     * @param unSubvolumeIdx List index of the subvolume.
     * @return Reference to the Subvolume object.
     */
    const detail::Subvolume & getSubvolume(UINTEGER unSubvolumeIdx) const
    {
      return const_cast<const detail::Subvolume &>(const_cast<DataModel *>(this)->getSubvolume(unSubvolumeIdx));
    };

    /**
     * Get a subvolume at a given position in the list.
     * 
     * @param unSubvolumeIdx List index of the subvolume.
     * @return Reference to the Subvolume object.
     */
    detail::Subvolume & getSubvolume(UINTEGER unSubvolumeIdx)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unSubvolumeIdx >= m_unSubvolumes)
        throw std::runtime_error("DataModel::getSubvolume() - invalid arguments.");
#endif
      return *(static_cast<detail::Subvolume *>(m_arSubvolumes[unSubvolumeIdx]));
    };

    /**
     * Get number of reaction wrappers.
     * 
     * @return Total number of Subvolume objects
     * associated with this DataModel.
     */
  inline UINTEGER getSubvolumesCount() const
    {
      return m_unSubvolumes;
    };

    /**
     * Print a string representation of the reaction network
     * 
     * @return The output stream
     */
    std::ostream & printReactionNetwork(std::ostream & os) const;

    //! Assignement operator
    DataModel& operator= (const DataModel&) = delete;

  ////////////////////////////////
  // Attributes
  protected:

    // Reactions
    //

    // Wrappers for reactions
    UINTEGER                        m_unReactionWrappers;  //!< Number of reaction wrappers
    detail::ReactionWrapper         *m_arReactionWrappers; //!< Array of reaction wrappers

    // Volumes
    //

    // Subvolumes
    UINTEGER                        m_unSubvolumes;   //!< Number of subvolumes
    detail::Subvolume               **m_arSubvolumes; //!< Array of subvolumes

    // Spatial dimensions
    BYTE                            m_uDims;     //!< Number of spatial dimensions
    UINTEGER                        *m_arunDims; //!< Array of dimension lengths

  ////////////////////////////////
  // Attributes
  public:
    // Species
    //

    // Reactions
    //

    //! Total propensity
    REAL                            dTotalPropensity;

    // Volumes
    //

    //! Information for sampling the next subvolume
    CompositionRejectionSamplerData crsdVolume;

  /////////////////////////////////////
  // Sampling variables
  public:
    //! Sampled reaction index
    INTEGER                         mu;

    //! Sampled volume
    UINTEGER                        nu;

    //! Sampled destination volume (for diffusion reactions)
    UINTEGER                        nu_D;

    //! Queued reactions (both D1 & D2 class reactions from (Cao et al, 2007))
    std::vector<DelayedReaction>    vQueuedReactions;
  };

}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_DATAMODEL_H_ */
