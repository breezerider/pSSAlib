/**
 * @file UpdateModule.h
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
 * Update module applies the changes to data structures due to reactions
 * using method-specific techniques
 */

#ifndef PSSALIB_UPDATE_UPDATEMODULE_H_
#define PSSALIB_UPDATE_UPDATEMODULE_H_

#include "../typedefs.h"

#ifdef PSSA_MODULE_LABEL
  #undef PSSA_MODULE_LABEL
#endif
#define PSSA_MODULE_LABEL pssalib::datamodel::SimulationInfo::eofModuleUpdate

// Forward declarations
namespace pssalib
{
namespace datamodel
{
  class SimulationInfo;
  
  namespace detail
  {
    class Subvolume;
    class Species;
    class SpeciesReference;
  }
} // close namespace datamodel

namespace update
{
  /**
   * @class UpdateModule
   * @brief This class serves as a base class for all update schemes 
   * offered by this library.
   * 
   * @details Provides methods to update method-specific data structures 
   * after a reaction has fired.
   * The respective methods of this class may be overloaded, offering a 
   * flexible way to extend its functionality.
   */
  class UpdateModule
  {
  /////////////////////////////////
  // Attributes
  protected:

    //! \internal Pointer to current reaction wrapper
    pssalib::datamodel::detail::ReactionWrapper
      * m_ptrReactionWrapper;

    // Pointers to subvolumes
    pssalib::datamodel::detail::Subvolume
      * m_ptrSubvolumeSrc, //!< pointer to source subvolume
      * m_ptrSubvolumeDst; //!< pointer to destination subvolume

    // Positions in the list
    UINTEGER m_sriBegin,     //!< index of first species reference for update
             m_sriEnd,       //!< index of last species reference for update
             m_sriReactants; //!< number of reactants in the list

  /////////////////////////////////
  // Constructors
  public:
    // Constructor
    UpdateModule();

    // Destructor
  virtual ~UpdateModule();

  /////////////////////////////////
  // Update module methods
  public:
    // Perform the update step
    bool doUpdate(pssalib::datamodel::SimulationInfo * ptrSimInfo);

  protected:
    // Schedule a delayed reaction
  virtual bool scheduleDelayed(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    //! Update volume structures
    bool updateVolumeStructures(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    // Update per species data structures after a chemical reaction
  virtual bool updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo) = 0;

    // Update per species data structures after a molecular diffusion event
  virtual bool updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo) = 0;
  };

}  } // close namespaces pssalib and update

#endif /* PSSALIB_UPDATE_UPDATEMODULE_H_ */
