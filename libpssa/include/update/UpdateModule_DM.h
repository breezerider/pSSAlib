/**
 * @file UpdateModule_DM.h
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
 * Update module definition for the Gillespie's Direct Method
 */

#ifndef PSSALIB_UPDATE_UPDATEMODULE_DM_H_
#define PSSALIB_UPDATE_UPDATEMODULE_DM_H_

#include "./UpdateModule.h"
#include "../../include/datamodel/DataModel_DM.h"

namespace pssalib
{
namespace datamodel
{
  class DataModel_DM;
  
  namespace detail
  {
    class Subvolume_DM;
  }
} // close namespace datamodel

namespace update
{
  /**
   * @class UpdateModule_DM
   * @brief Update the datastructures for the Gillespie's Direct Method 
   * with changes due to the fired reaction.
   *
   * @copydetails UpdateModule
   */
  class UpdateModule_DM : public UpdateModule
  {
  ////////////////////////////////
  // Constructors
  public:
    // Constructor
    UpdateModule_DM();

    // Destructor
virtual ~UpdateModule_DM();

  ////////////////////////////////
  // Update module methods
  protected:
    //! Update per species data structures after a chemical reaction
virtual bool updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    //! Update per species data structures after a molecular diffusion event
virtual bool updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    //! Update data structures
    bool updateSpeciesStructures(pssalib::datamodel::DataModel_DM * ptrDMData,
                                 pssalib::datamodel::detail::Subvolume_DM & DMSubVol);
  };

}  } // close namespaces pssalib and update

#endif /* PSSALIB_UPDATE_UPDATEMODULE_DM_H_ */
