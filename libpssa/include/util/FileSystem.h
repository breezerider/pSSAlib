/**
 * @file FileSystem.h
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
 * Auxiliary functions for handling requests to the file system
 */

#ifndef PSSALIB_UTIL_FILESYSTEM_H_
#define PSSALIB_UTIL_FILESYSTEM_H_

#include "../typedefs.h"
#include "../stdheaders.h"

namespace pssalib
{
namespace util
{
  //! \internal A cross-platform wrapper for mkdir function
  bool makeDir(std::vector< STRING > &, STRING &, bool dryRun = false);

  //! \internal A cross-platform wrapper for mkdir function
  bool makeDir(STRING &, bool dryRun = false);

  //! \internal A cross-platform directory path validator
  bool checkPath(const STRING &);

  //! \internal A cross-platform file path constructor
  void makeFilePath(const STRING &, const STRING &, STRING &);

} } // close namespaces util and pssalib

#endif /* PSSALIB_UTIL_FILESYSTEM_H_ */
