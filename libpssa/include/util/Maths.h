/**
 * @file Maths.h
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
 * Common header file with definitions for maths functions.
 */

#ifndef PSSALIB_UTIL_MATHS_H_
#define PSSALIB_UTIL_MATHS_H_

// Log base 2
#ifndef LOG2
# if defined(_WIN32)
#  define LOG2(_X) (log(_X) * M_LOG2E/M_LOG10E)
# else // hope this works on all other platforms
#  define LOG2(_X) log2(_X)
# endif
#endif

namespace pssalib
{
namespace maths
{

  /**
   * An optimized version of floor(log2(x)) for integer arguments
   * Based on http://aggregate.org/MAGIC/
   * 
   * @param x argument
   * @return floor(log2(x))
   */
  inline unsigned int floor_log2(unsigned int x)
  {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    x >>= 1;

    x -= ((x >> 1) & 0x55555555);
    x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
    x = (((x >> 4) + x) & 0x0f0f0f0f);
    x += (x >> 8);
    x += (x >> 16);
    return(x & 0x0000003f);
  }

} } /* close namespaces maths and pssalib */

#endif /* PSSALIB_UTIL_MATHS_H_ */
