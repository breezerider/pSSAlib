/**
 * @file SamplingModule.h
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
 * Sampling module acquires random samples using method-specific techniques
 */

#ifndef PSSALIB_SAMPLING_SAMPLINGMODULE_H_
#define PSSALIB_SAMPLING_SAMPLINGMODULE_H_

#include "../typedefs.h"
#include "./CompositionRejectionSampler.h"

#ifdef PSSA_MODULE_LABEL
  #undef PSSA_MODULE_LABEL
#endif
#define PSSA_MODULE_LABEL pssalib::datamodel::SimulationInfo::eofModuleSampling

// Forward declarations
namespace pssalib
{
namespace datamodel
{
  class SimulationInfo;
} // close namespace datamodel

namespace sampling
{
  /**
   * @class SamplingModule
   * @brief This class serves as a base class for all sampling schemes 
   * offered by this library.
   * 
   * @details Provides methods to sample at random in method-specific manner.
   * The respective methods of this class may be overloaded, offering a 
   * flexible way to extend its functionality.
   */
  class SamplingModule
  {
  ////////////////////////////////
  // Attributes
  protected:
    //! Pseudo-Random numbers generator from GSL
    gsl_rng * m_ptrRNG;
    CompositionRejectionSampler crVolumeSampler;

  ////////////////////////////////
  // Constructors
  public:
    // Default Constructor
    SamplingModule();

    // Destructor
    virtual ~SamplingModule();

  //////////////////////////////
  // Methods
  protected:
    // Sample next reaction time
    bool sampleTime(pssalib::datamodel::SimulationInfo* ptrSimInfo);

    // Sample subvolume
    bool sampleVolume(pssalib::datamodel::SimulationInfo* ptrSimInfo);

    //! Sample next reaction index
    virtual bool sampleReaction(pssalib::datamodel::SimulationInfo* ptrSimInfo) = 0;

  public:
    // Set the seed of the random number generator
    void set_rng_seed(UINTEGER seed);

    // Get next sample
    bool getSample(pssalib::datamodel::SimulationInfo* ptrSimInfo);
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_SAMPLINGMODULE_H_ */
