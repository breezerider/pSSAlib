/**
 * @file UpdateModule_SPDM.h
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
 * Update module definition for the Sorting Partial Propensity Direct Method
 */

#ifndef PSSALIB_UPDATE_UPDATEMODULE_SPDM_H_
#define PSSALIB_UPDATE_UPDATEMODULE_SPDM_H_

#include "./UpdateModule_PDM.h"

namespace pssalib
{
namespace update
{
  /**
   * @class UpdateModule_SPDM
   * @brief Update the datastructures for the Sorting Partial 
   * Propensity Direct Method with changes due to the fired reaction.
   *
   * @copydetails UpdateModule
   */
  class UpdateModule_SPDM : public UpdateModule_PDM
  {
  ////////////////////////////////
  // Constructors
  public:
    // Default Constructor
    UpdateModule_SPDM();

    // Destructor
virtual ~UpdateModule_SPDM();

  ////////////////////////////////
  // Update module methods
  protected:
    //! Update per species data structures after a chemical reaction
virtual bool updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo);
  };

}  } // close namespaces pssalib and update

#endif /* PSSALIB_UPDATE_UPDATEMODULE_SPDM_H_ */
