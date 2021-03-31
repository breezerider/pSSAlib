/**
 * @file SamplingModule_DM.cpp
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
 * Sampling module implementation for the Gillespie's Direct Method
 */

#include "../../include/datamodel/DataModel_DM.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/sampling/SamplingModule_DM.h"

namespace pssalib
{
namespace sampling
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  SamplingModule_DM::SamplingModule_DM()
  {
    // Do nothing
  }

  //! Destructor
  SamplingModule_DM::~SamplingModule_DM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Sample reaction index
  bool SamplingModule_DM::sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_DM * ptrDMData = static_cast<pssalib::datamodel::DataModel_DM*>(
      ptrSimInfo->getDataModel());
    pssalib::datamodel::detail::Subvolume_DM & DMSubVol =
      ptrDMData->getSubvolume(ptrDMData->nu);

    UINTEGER mu,
             M = ptrDMData->getReactionWrappersCount();

    // Sample reaction
    REAL temp1 = gsl_rng_uniform_pos (m_ptrRNG) * DMSubVol.dTotalPropensity;
    REAL temp2 = 0.0;
    for(mu = 0; mu < M; ++mu)
    {
      temp2 += DMSubVol.propensity(mu);

      if(temp1 <= temp2)
        break;
    }

    ptrDMData->mu = mu;

    return true;
  }

}  } // close namespaces pssalib and sampling
