/**
 * @file SamplingModule_SPDM.h
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
 * Sampling module definition for the Sorting Partial Propensity 
 * Direct Method
 */

#ifndef PSSALIB_SAMPLING_SAMPLINGMODULE_SPDM_H_
#define PSSALIB_SAMPLING_SAMPLINGMODULE_SPDM_H_

#include "./SamplingModule_PDM.h"

namespace pssalib
{
namespace sampling
{
  /**
   * \class SamplingModule_SPDM
   * \brief Provide random samples using the Sorting Partial Propensity 
   * Direct Method.
   * 
   * \copydetails SamplingModule
   */
  class SamplingModule_SPDM : public SamplingModule_PDM
  {
  /////////////////////////////////
  // Constructors
  public:
    // Constructor
    SamplingModule_SPDM();
    // Destructor
    virtual ~SamplingModule_SPDM();

  ////////////////////////////////
  // Methods
  protected:
    // Sample reaction index
    virtual bool sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_SAMPLINGMODULE_SPDM_H_ */
