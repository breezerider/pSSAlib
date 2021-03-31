/**
 * @file SamplingModule_PDM.h
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
 * Sampling module definition for the Partial Propensity 
 * Direct Method
 */

#ifndef PSSALIB_SAMPLING_SAMPLINGMODULE_PDM_H_
#define PSSALIB_SAMPLING_SAMPLINGMODULE_PDM_H_

#include "./SamplingModule.h"

namespace pssalib
{
namespace sampling
{
  /**
   * @class SamplingModule_PDM
   * @brief Provide random samples using the Partial Propensity Direct Method.
   * 
   * @copydetails SamplingModule
   */
  class SamplingModule_PDM : public SamplingModule
  {
  /////////////////////////////////
  // Constructors
  public:
    // Constructor
    SamplingModule_PDM();
    // Destructor
  virtual ~SamplingModule_PDM();

  ////////////////////////////////
  // Methods
  protected:
    // Sample reaction index
    virtual bool sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_SAMPLINGMODULE_PDM_H_ */
