/**
 * @file InplaceMemory.h
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
 * Auxiliary functions for allocating objects in pre-allocated memory buffers
 */

#ifndef PSSALIB_UTIL_INPLACE_MEMORY_H_
#define PSSALIB_UTIL_INPLACE_MEMORY_H_

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
  template< class T >
  T * inplace_alloc(const UINTEGER n, const T & t)
  {
    if(0 == n)
      return NULL;

    unsigned char * buffer = new unsigned char[n*sizeof(T)];

    for(UINTEGER i = 0; i < n; i++)
      new (buffer + i*sizeof(T)) T(t);

    return (T *) buffer;
  }

  template< class T >
  T * inplace_alloc(const UINTEGER n)
  {
    if(0 == n)
      return NULL;

    unsigned char * buffer = new unsigned char[n*sizeof(T)];

    return (T *) buffer;
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
  template< class T >
  void inplace_free(UINTEGER & n, T ** t)
  {
    if(0 == n)
      return;

    T * temp = *t;
    for(UINTEGER i = 0; i < n; i++)
      (temp++)->T::~T();

    delete [] ((unsigned char *)(*t));
    (*t) = NULL;
    n = 0;
  }

} } // close namespaces util and pssalib

#endif /* PSSALIB_UTIL_INPLACE_MEMORY_H_ */
