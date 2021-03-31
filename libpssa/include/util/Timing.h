/**
 * @file Timing.h
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
 * Auxiliary functions for producing time points
 */

#ifndef PSSALIB_UTIL_TIMING_H_
#define PSSALIB_UTIL_TIMING_H_

#include "../typedefs.h"

namespace pssalib
{
namespace timing
{
  /**
   * Return number of time points within an interval
   * 
   * @param tb First check point in simulation output
   * @param te Last check point in simulation output
   * @param dt Distance between check points
   * @return Number of time check point
   */
  inline UINTEGER getNumTimePoints(REAL tb, REAL te, REAL dt)
  {
    if((tb >= 0.0)&&(tb < te)&&(dt > 0.0))
    {
// std::cout << "tb = " << tb << "; te = " << te << "; dt = " << dt << "\n";
//       UINTEGER unTPS = 2;
//       REAL t = te - tb;
//       while(t > dt)
//       {
//         t -= dt;
//         unTPS++;
// std::cout << "t = " << t << "; unTPS = " << unTPS << "\n";
//       }
//       return unTPS;
      const REAL nom = (te - tb);
      UINTEGER unTPS = std::floor(nom/dt);
      const REAL t = tb + REAL(unTPS)*dt;

      if(t == te)
        return unTPS+1;
      else
        return unTPS+2;
    }
    return 0;
  }

  /**
   * Return the check-point that is greater or equal
   * to a given moment in simulation time
   * 
   * @param tb First check point in simulation output
   * @param te Last check point in simulation output
   * @param dt Distance between check points
   * @param t  Any time point in simulation time
   * @return A time check point
   */
  inline REAL getAdjTimePointHi(REAL tb, REAL te, REAL dt, REAL t)
  {
    if((tb >= 0.0)&&(tb < te)&&(dt > 0.0)&&(t >= tb))
    {
      t = tb + ((REAL)std::ceil((t - tb) / dt)) * dt;

      return std::max(std::min(t, te), tb);
    }
    return std::numeric_limits<REAL>::infinity();
  }

  /**
   * Return the check-point that is less or equal
   * to a given moment in simulation time
   * 
   * @param tb First check point in simulation output
   * @param te Last check point in simulation output
   * @param dt Distance between check points
   * @param t  Any time point in simulation time
   * @return A time check point
   */
  inline REAL getAdjTimePointLo(REAL tb, REAL te, REAL dt, REAL t)
  {
    if((tb >= 0.0)&&(tb < te)&&(dt > 0.0)&&(t >= tb))
    {
      t = tb + ((REAL)std::floor((t - tb) / dt)) * dt;

      return std::max(std::min(t, te), tb);
    }
    return REAL(-1.0) * std::numeric_limits<REAL>::infinity();
  }

} } // close timing & pssalib namespaces

#endif /* PSSALIB_UTIL_TIMING_H_ */
