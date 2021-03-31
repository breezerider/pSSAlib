/**
 * @file FileSystem.cpp
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
 * Implementation of auxiliary functions for handling requests to the file system
 */

#include "../../include/util/FileSystem.h"
#include "../../include/util/MPIWrapper.h"

#include <boost/algorithm/string.hpp>

namespace pssalib
{
namespace util
{
  static const STRING::value_type path_separator =
#ifdef _WIN32
    '\\';
#else /* *NIX, MACOS */
    '/';
#endif
  
  //! @internal A cross-platform wrapper for mkdir function
  bool makeDir(std::vector< STRING > & arPath, STRING & outPath, bool dryRun)
  {
    // Result
    bool bResult = false;

    // prepare the output variable
    outPath.clear();

    // validate arguments
    if(!arPath.empty())
    {
      // loop vars
      for(std::vector< STRING >::iterator it = arPath.begin(); it != arPath.end(); ++it) 
      {
        STRING tmp(*it);
        boost::algorithm::trim(tmp);
        if(tmp.empty())
          continue;
        if(path_separator != *(tmp.rbegin()))
          tmp.append(1, path_separator);

        outPath += tmp;
      }

      if(PSSALIB_MPI_IS_MASTER)
        bResult = makeDir(outPath, dryRun);
#ifdef HAVE_MPI
      bResult = getMPIWrapperInstance().broadcast(&bResult, sizeof(bool), 0);
#endif
    }

    return bResult;
  }

  //! @internal A cross-platform wrapper for mkdir function
  bool makeDir(STRING & path, bool dryRun)
  {
    // Result
    bool bResult = false;

    // validate arguments
    boost::algorithm::trim(path);
    if (!path.empty())
    {
      bResult = true;

      // loop vars
      STRING::size_type pos = 0;

      // append a path separator at the end
      if (path_separator != *(path.rbegin()))
        path.append(1, path_separator);

      STRING::size_type len = path.length();
      while(pos < (len-1))
      {
        pos = path.find(path_separator, pos + 1);
        if (pos == STRING::npos)
        {
          bResult = false;
          break;
        }

        if(PSSALIB_MPI_IS_MASTER&&!dryRun)
        {
          if(
#ifdef _WIN32
            (_tmkdir(path.substr(0, pos).c_str()) != 0)
#else /* *NIX, MACOS */
            (mkdir(path.substr(0, pos).c_str(),
              S_IRWXU | S_IRWXG | S_IRWXO) != 0)
#endif
            &&(errno != EEXIST))
          {
            bResult = false;
            break;
          }
        }
      }
    }

    return bResult;
  }

  //! @internal A cross-platform directory path validator
  bool checkPath(const STRING &outPath)
  {
#ifdef _WIN32
    return ((TRUE == PathFileExists(outPath.c_str()))&&
            (TRUE == PathIsDirectory(outPath.c_str())));
#else /* *NIX, MACOS */
    struct stat sb;
    return ((0 == stat(outPath.c_str(), &sb))&&
            S_ISDIR(sb.st_mode));
#endif
  }

  //! @internal A cross-platform file path constructor
  void makeFilePath(const STRING & basePath, const STRING & fileName, STRING & outPath)
  {
    // prepare the output variable
    outPath.assign(basePath);

    // append a path separator at the end
    if(path_separator != *(outPath.rbegin()))
      outPath.append(1, path_separator);
    outPath.append(fileName);
  }
} } // close namespaces util and pssalib
