/**
 * @file UpdateModule_SPDM.cpp
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
 * Update module implementation for the Sorting Partial Propensity
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/update/UpdateModule_SPDM.h"

namespace pssalib
{
namespace update
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  UpdateModule_SPDM::UpdateModule_SPDM()
  {
    // Do nothing
  }

  //! Destructor
  UpdateModule_SPDM::~UpdateModule_SPDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods
  
  bool UpdateModule_SPDM::updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    if(!UpdateModule_PDM::updateSpeciesStructuresReaction(ptrSimInfo))
      return false;

    pssalib::datamodel::DataModel_SPDM * ptrSPDMData = 
      static_cast<pssalib::datamodel::DataModel_SPDM * >
        (ptrSimInfo->getDataModel());

    if(ptrSPDMData->rowIndex > 0)
      // Swap with preceding
      static_cast<pssalib::datamodel::detail::Subvolume_SPDM *>(m_ptrSubvolumeSrc)->moveRowUp(ptrSPDMData->rowIndex);
    if(ptrSPDMData->colIndex > 0)
      // Swap with preceding
      static_cast<pssalib::datamodel::detail::Subvolume_SPDM *>(m_ptrSubvolumeSrc)->moveColLeft(ptrSPDMData->rowIndex, ptrSPDMData->colIndex);

    return true;
  }

}  } // close namespaces pssalib and update
