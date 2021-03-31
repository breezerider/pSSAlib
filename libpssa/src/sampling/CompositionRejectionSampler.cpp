/**
 * @file CompositionRejectionSampler.cpp
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
 * Composition-Rejection Sampler Implementation
 */

#include "../../include/sampling/CompositionRejectionSampler.h"
#include "../../include/datamodel/CompositionRejectionSamplerData.h"

namespace pssalib
{
namespace sampling
{
  bool CompositionRejectionSampler::Sample(const pssalib::datamodel::CompositionRejectionSamplerData * ptrData, 
                                           gsl_rng* ptrRNG, const REAL scale, UINTEGER& outI, REAL& outR)
  {
    for(UINTEGER k = 0; k < PSSA_CR_MAX_ITER; ++k)
    {
      REAL r = gsl_rng_uniform (ptrRNG) * scale;

      pssalib::datamodel::PSSACR_Bins::CONST_PAIR_BINS_ITER itB;
      int c = 0;
      ptrData->bins.getBins(itB);

      // Linear search step to find the bin
      REAL temp = 0.0;
      for(; itB.first != itB.second; ++itB.first)
      {
        c++;

        temp += itB.first->second.dBinSum;
        if(r < temp)
          break;
      }

      // FIXME no particular order is guaranteed when traversing boost::unordered_map
      // When r =~ scale, In some cases the sum can fail for small value and reach the end without
      // find a bin. In this case take the latest one bin
      if (itB.first == itB.second)
      {
        // It seem that boost does not support --:
        // so reiterate c - 1

        ptrData->bins.getBins(itB);

        for(int s = 0; s <  c - 1 ; s++) ++itB.first;
      }

      if (itB.first->first <= 30)
        temp = ptrData->minValue * (1 << itB.first->first);
      else
        temp = ldexp(ptrData->minValue, itB.first->first);
      const pssalib::datamodel::PSSACR_Bin *pBin = &(itB.first->second);
      UINTEGER unBins = pBin->size();

      // Rejection step to sample within the bin
      while (true)
      {
        UINTEGER sI = gsl_rng_uniform_int(ptrRNG, unBins);
        r = gsl_rng_uniform (ptrRNG);
        r *= temp;

        sI = pBin->get_at(sI);

        if(r < ptrData->bins.getValue(sI)) {
          outI = sI;
          outR = r;
          return true;
        }
      }
    }

    // All bins are guaranteed to be at least 50% full, so the probability of 
    // not finding a sample is <= 2^(- PSSA_CR_MAX_ITER) (~8e-31).
    return false;
  }

}  } // close namespaces pssalib and sampling
