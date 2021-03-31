/**
 * @file Indexing.h
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
 * Auxiliary functions for computing
 */

#ifndef PSSALIB_UTIL_INDEXING_H_
#define PSSALIB_UTIL_INDEXING_H_

#include "../typedefs.h"
#include "../stdheaders.h"

namespace pssalib
{
namespace util
{
  /**
   * Converts linear index of an array entry
   * to the corresponding subscript indexes.
   * 
   * @param uDims Number of array dimensions.
   * @param arunDims Array of dimension lengths.
   * @param unIdx Linear index of the array entry.
   * @param arunSub A pre-allocated output array
   * for the corresponding subscript indexes.
   */
  inline void ind2sub(BYTE uDims, const UINTEGER * arunDims, UINTEGER unIdx, UINTEGER * arunSub)
  {     
    UINTEGER temp1 = unIdx, temp2 = unIdx;
    for(BYTE di = 0; di < uDims - 1; ++di)
    {
        temp2 = temp1 / arunDims[di];
        arunSub[di] = temp1 - temp2 * arunDims[di];
        temp1 = temp2;
    }
    arunSub[uDims-1] = temp2;
  }

  /**
   * Converts subscript indexes of an array
   * entry to the corresponding linear index.
   * 
   * @param uDims Number of array dimensions.
   * @param arunDims Array of dimension lengths.
   * @param arunSub Array of subscript indexes.
   * @param unIdx Linear index of the array entry.
   */
  inline void sub2ind(BYTE uDims, const UINTEGER * arunDims, UINTEGER * arunSub, UINTEGER & unIdx)
  {     
    UINTEGER temp1 = 1, temp2 = 0;
    for(BYTE di = 0; di < uDims; ++di)
    {
        temp2 += arunSub[di] * temp1;
        temp1 *= arunDims[di];
    }
    unIdx = temp2;
  }

} } // close namespaces util and pssalib

#endif /* PSSALIB_UTIL_INDEXING_H_ */
