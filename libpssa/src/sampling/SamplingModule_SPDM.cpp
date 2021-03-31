/**
 * @file SamplingModule_PDM.cpp
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
 * Sampling module implementation for the Sorting Partial Propensity 
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/sampling/SamplingModule_SPDM.h"

namespace pssalib
{
namespace sampling
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  SamplingModule_SPDM::SamplingModule_SPDM()
  {
    // Do nothing
  }

  //! Destructor
  SamplingModule_SPDM::~SamplingModule_SPDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Sampling module methods

  //! Sample next reaction index using SpDM
  bool SamplingModule_SPDM::sampleReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo)
#define PSSALIB_INTERNAL_SPDM_MODULE
#include "./inc/SamplingModule_S_PDM.inc"
#undef PSSALIB_INTERNAL_SPDM_MODULE

}  } // close namespaces pssalib and sampling
