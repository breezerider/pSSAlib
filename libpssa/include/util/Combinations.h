/**
 * @file Combinations.h
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
 * Auxiliary functions for computing numbers of reaction combinations
 */

#ifndef PSSALIB_UTIL_COMBINATIONS_H_
#define PSSALIB_UTIL_COMBINATIONS_H_

#include "../typedefs.h"
#include "../stdheaders.h"

namespace pssalib
{
namespace util
{
  /**
   * Get the number of possible reaction combinations.
   * 
   * @param n number of molecules of reactant species.
   * @param m reactant stoichiometry.
   * @return number of possible reaction combinations.
   */
  inline REAL getPartialCombinationsHeteroreactions(ULINTEGER n, ULINTEGER m)
  {
    if(n < m)
      return 0.0;
    // most frequent cases
    else if(1 == m)
      return REAL(n);
    else if(2 == m)
      return REAL((n * (n - 1)) / 2);
    else if(0 == m)
      return 1.0;

    ULINTEGER nom = --n, den = m--;
    for(; m > 1; --m)
    {
      nom *= --n; den *= m;
    }

    return ((REAL)nom*n)/((REAL)den);
  }

  /**
   * Get the number of possible reaction combinations
   * (designed for partial-propensity methods when
   * considering uni-molecular and self-dependent reactions).
   * 
   * @param n number of molecules of reactant species.
   * @param m reactant stoichiometry.
   * @return number of possible reaction combinations.
   */
  inline REAL getPartialCombinationsHomoreactions(ULINTEGER n, ULINTEGER m)
  {
    // uni-molecular case
    if((1 == m)||(0 == m))
      return 1.0;
    else
    {
      if(n < m)
        return 0.0;
      // most frequent case
      else if(2 == m)
        return REAL(--n) / REAL(2.0);
      else
      {
        ULINTEGER nom = --n, den = m--;
        for(; m > 1; --m)
        {
          nom *= --n; den *= m;
        }

        return ((REAL)nom)/((REAL)den);
      }
    }
  }

} } // close namespaces util and pssalib

#endif /* PSSALIB_UTIL_COMBINATIONS_H_ */
