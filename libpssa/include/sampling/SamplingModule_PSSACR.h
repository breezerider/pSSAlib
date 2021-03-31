/**
 * @file SamplingModule_PSSACR.h
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
 * Sampling module definition for the Partial Propensity SSA with
 * Composition-Rejection Sampling
 */

#ifndef PSSALIB_SAMPLING_SAMPLINGMODULE_PSSACR_H_
#define PSSALIB_SAMPLING_SAMPLINGMODULE_PSSACR_H_

#include "./SamplingModule_PDM.h"
#include "./CompositionRejectionSampler.h"

namespace pssalib
{
namespace sampling
{
  /**
   * @class SamplingModule_PSSACR
   * @brief Provide random samples using the Partial Propensity Direct Method 
   * with Compositio-Rejection sampling.
   * 
   * @copydetails SamplingModule
   */
  class SamplingModule_PSSACR : public SamplingModule_PDM
  {
  ////////////////////////////////
  // Constructors
  public:
    // Default Constructor
    SamplingModule_PSSACR();
    // Destructor
    virtual	~SamplingModule_PSSACR();

  //////////////////////////////
  // Methods
  protected:
    // Sample reaction index
    virtual bool sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);

  ////////////////////////////////
  // Attributes
  protected:
    CompositionRejectionSampler crSampler;
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_SAMPLINGMODULE_PSSACR_H_ */
