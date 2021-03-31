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
 * Sampling module implementation for the Partial-Propensity 
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/sampling/SamplingModule_PDM.h"

namespace pssalib
{
namespace sampling
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  SamplingModule_PDM::SamplingModule_PDM()
  {
    // Do nothing
  }

  //! Destructor
  SamplingModule_PDM::~SamplingModule_PDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Sampling module methods

  //! Sample next reaction index using pDM
  bool SamplingModule_PDM::sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
#define PSSALIB_INTERNAL_PDM_MODULE
#include "./inc/SamplingModule_S_PDM.inc"
#undef PSSALIB_INTERNAL_PDM_MODULE

}  } // close namespaces pssalib and sampling
