/**
 * @file SamplingModule_DM.h
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
 * Sampling module definition for the Gillespie's Direct Method
 */

#ifndef PSSALIB_SAMPLING_SAMPLINGMODULE_DM_H_
#define PSSALIB_SAMPLING_SAMPLINGMODULE_DM_H_

#include "./SamplingModule.h"

namespace pssalib
{
namespace sampling
{
  /**
   * @class SamplingModule_DM
   * @brief Provide random samples using the Gillespie's Direct Method.
   * 
   * @copydetails SamplingModule
   */
  class SamplingModule_DM : public SamplingModule
  {
  /////////////////////////////////////
  // Constructors
  public:
    // Default Constructor
    SamplingModule_DM();
    // Destructor
  virtual ~SamplingModule_DM();

  //////////////////////////////
  // Methods
  protected:
    // Sample next reaction index
    virtual bool sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo);
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_SAMPLINGMODULE_DM_H_ */
