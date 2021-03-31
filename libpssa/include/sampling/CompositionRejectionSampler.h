/**
 * @file CompositionRejectionSampler.h
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
 * A stand-alone Sampler implementation using the Composition-Rejection method
 */

#ifndef PSSALIB_SAMPLING_COMPOSITION_REJECTION_SAMPLER_H_
#define PSSALIB_SAMPLING_COMPOSITION_REJECTION_SAMPLER_H_

#include "../typedefs.h"

namespace pssalib
{
// Forward declarations
namespace datamodel
{
  class CompositionRejectionSamplerData;
}

namespace sampling
{
  class CompositionRejectionSampler
  {
    public:
      bool Sample(const pssalib::datamodel::CompositionRejectionSamplerData* ptrData, 
                  gsl_rng* ptrRNG, const REAL scale, UINTEGER& outI, REAL& outR);
  };

}  } // close namespaces pssalib and sampling

#endif /* PSSALIB_SAMPLING_COMPOSITION_REJECTION_SAMPLER_H_ */

