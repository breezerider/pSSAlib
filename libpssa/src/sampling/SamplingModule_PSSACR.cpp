/**
 * @file SamplingModule_PSSACR.cpp
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
 * Sampling module implementation for the Partial Propensity SSA with
 * Composition-Rejection Sampling (Ramaswamy, 2010)
 */

#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/datamodel/PSSACR_Bins.h"
#include "../../include/sampling/SamplingModule_PSSACR.h"

namespace pssalib
{
namespace sampling
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  SamplingModule_PSSACR::SamplingModule_PSSACR()
  {
    // Do nothing
  }

  //! Destructor
  SamplingModule_PSSACR::~SamplingModule_PSSACR()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Sampling module methods

  //! Sample next reaction index using pSSA-CR
  bool SamplingModule_PSSACR::sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_PSSACR * ptrPSRDCRData =
      static_cast<pssalib::datamodel::DataModel_PSSACR *>(ptrSimInfo->getDataModel());
  pssalib::datamodel::detail::Subvolume_PSSACR & PSSACRSubVol = ptrPSRDCRData->getSubvolume(ptrPSRDCRData->nu);

    REAL r = 0.0;
    UINTEGER sI = 0;
    bool success = crSampler.Sample(&PSSACRSubVol.crsdSigma, m_ptrRNG, PSSACRSubVol.dTotalPropensity, sI, r);

    if (success) {
      UINTEGER sJ = 0;
      success = crSampler.Sample(&PSSACRSubVol.crsdPi(sI), m_ptrRNG, PSSACRSubVol.lambda(sI), sJ, r);

      if (success) {
        ptrPSRDCRData->mu = ptrPSRDCRData->aruL(sI, sJ)->getSerialNumber();
      }
    }

    if (!success) {
      PSSA_ERROR(ptrSimInfo, << "sampling did not converge in given "
        "number of iterations." << std::endl);
    }
    return success;
  }

}  } // close namespaces pssalib and sampling
