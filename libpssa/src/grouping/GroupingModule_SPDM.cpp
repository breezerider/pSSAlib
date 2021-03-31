/**
 * @file GroupingModule_SPDM.cpp
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
 * Grouping implementation for the Sorting Partial Propensity
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/grouping/GroupingModule_SPDM.h"

namespace pssalib
{
namespace grouping
{
  ////////////////////////////////
  // Constructors

  //! Constructor
  GroupingModule_SPDM::GroupingModule_SPDM()
  {
    // Do nothing
  }

  //! Copy constructor
  GroupingModule_SPDM::GroupingModule_SPDM(GroupingModule &g) 
    : GroupingModule_PDM(g)
  {
    // Do nothing
  }

  //! Destructor
  GroupingModule_SPDM::~GroupingModule_SPDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Initialize data structures (called before each trial)
  bool GroupingModule_SPDM::initialize(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    // Call the baseclass method
    if(!GroupingModule_PDM::initialize(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_SPDM * ptrSPDMData = 
      static_cast<pssalib::datamodel::DataModel_SPDM *>
        (ptrSimInfo->getDataModel());

    // Reset indexing
    for(UINTEGER svi = 0; svi < ptrSPDMData->getSubvolumesCount(); ++svi)
      ptrSPDMData->getSubvolume(svi).resetIndexing();
    
    ptrSPDMData->rowIndex = ptrSPDMData->colIndex = 0;

    return true;
  }

}  } // close namespaces pssalib and grouping
