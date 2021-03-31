/**
 * @file UpdateModule_PDM.h
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
 * Update module definition for the Partial-Propensity Direct Method
 */

#ifndef PSSALIB_UPDATE_UPDATEMODULE_PDM_H_
#define PSSALIB_UPDATE_UPDATEMODULE_PDM_H_

#include "./UpdateModule.h"

namespace pssalib
{
namespace datamodel
{
  class DataModel_PDM;

  namespace detail
  {
    class Subvolume_PDM;
  }
} // close namespace datamodel

namespace update
{
  /**
   * @class UpdateModule_PDM
   * @brief Update the data structures for the Partial Propensity Direct 
   * Method with changes due to the fired reaction.
   *
   * @copydetails UpdateModule
   */
  class UpdateModule_PDM : public UpdateModule
  {
  ////////////////////////////////
  // Constructors
  public:
    // Default Constructor
    UpdateModule_PDM();

    // Destructor
    virtual ~UpdateModule_PDM();

  ////////////////////////////////
  // Update module methods
  protected:
    // Update per species data structures after a chemical reaction
virtual bool updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    // Update per species data structures after a molecular diffusion event
virtual bool updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    // Update data structures
    bool updateSpeciesStructures(pssalib::datamodel::SimulationInfo * ptrSimInfo,
                                 pssalib::datamodel::DataModel_PDM * ptrPDMData,
                                 pssalib::datamodel::detail::Subvolume_PDM & PDMSubVol,
                                 UINTEGER index, INTEGER stoichiometry);
  };

}  } // close namespaces pssalib and grouping

#endif /* PSSALIB_UPDATE_UPDATEMODULE_PDM_H_ */
