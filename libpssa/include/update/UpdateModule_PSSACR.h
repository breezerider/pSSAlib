/**
 * @file UpdateModule_PSSACR.h
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
 * Update module definition for the Partial-Propensity SSA 
 * with Composition-Rejection sampling
 */

#ifndef PSSALIB_UPDATE_UPDATEMODULE_PSSACR_H_
#define PSSALIB_UPDATE_UPDATEMODULE_PSSACR_H_

#include "./UpdateModule_PDM.h"

namespace pssalib
{
namespace datamodel
{
  class DataModel_PSSACR;
  
  namespace detail
  {
    class Subvolume_PSSACR;
  }
} // close namespace datamodel

namespace update
{
  /**
   * @class UpdateModule_PSSACR
   * @brief Update the datastructures for the Partial Propensity Direct 
   * Method with Composition-Rejection sampling with changes due to the 
   * fired reaction.
   *
   * @copydetails UpdateModule
   */
  class UpdateModule_PSSACR : public UpdateModule_PDM
  {
  ////////////////////////////////
  // Constructors
  public:
    // Constructor
    UpdateModule_PSSACR();

    // Destructor
virtual ~UpdateModule_PSSACR();

  ////////////////////////////////
  // Methods
  protected:
    //! Update per species data structures after a chemical reaction
virtual bool updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    //! Update per species data structures after a molecular diffusion event
virtual bool updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    //! Update per species data structures
    void updateSpeciesStructures(pssalib::datamodel::DataModel_PSSACR * ptrPSSACRData,
                                 pssalib::datamodel::detail::Subvolume_PSSACR & subVol,
                                 const UINTEGER index);
  };

}  } // close namespaces pssalib and update

#endif /* PSSALIB_UPDATE_UPDATEMODULE_PSSACR_H_ */
