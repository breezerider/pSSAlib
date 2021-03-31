/**
 * @file typedefs.h
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
 * Common header file with type definitions.
 */

#include "./stdheaders.h"

#ifndef PSSALIB_TYPEDEFS_H_
#define PSSALIB_TYPEDEFS_H_

///////////////////////////////////////////////////////////////
// Common output macros

#ifdef DEBUG
#  define PSSALIB_INTERNAL_LOG(si, type, header, args) \
          while((si)->isLoggingOn(type)) { \
            (si)->report() << BOOST_CURRENT_FUNCTION << " " << header args; break; \
          };
#else
#  define PSSALIB_INTERNAL_LOG(si, type, header, args) \
          while((si)->isLoggingOn(type)) { \
            (si)->report() << header args; break; \
          };
#endif

#define PSSA_TRACE(si, args) \
          PSSALIB_INTERNAL_LOG(si, pssalib::datamodel::SimulationInfo::ofTrace | PSSA_MODULE_LABEL, "(TRACE) : ", args)

#define PSSA_INFO(si, args) \
          PSSALIB_INTERNAL_LOG(si, pssalib::datamodel::SimulationInfo::ofInfo, "(INFO) : ", args)

#define PSSA_WARNING(si, args) \
          PSSALIB_INTERNAL_LOG(si, pssalib::datamodel::SimulationInfo::ofWarning, "(WARNING) : ", args)

#define PSSA_ERROR(si, args) \
          PSSALIB_INTERNAL_LOG(si, pssalib::datamodel::SimulationInfo::ofError, "(ERROR) : ", args)

///////////////////////////////////////////////////////////////

#define PSSA_CR_MAX_ITER 100
#define PSSA_MODULE_LABEL pssalib::datamodel::SimulationInfo::ofNone

#define PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER ","
#define PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER "\t"

///////////////////////////////////////////////////////////////
// Define some basic types
#ifndef BYTE_DEFINED
  #undef BYTE
  typedef unsigned char  BYTE;
  #define BYTE_DEFINED
#endif
#ifndef SHORT_DEFINED
  #undef SHORT
  typedef short          SHORT;
  #define SHORT_DEFINED
#endif
#ifndef USHORT_DEFINED
  #undef USHORT
  typedef unsigned short USHORT;
  #define USHORT_DEFINED
#endif
#ifndef INTEGER_DEFINED
  #undef INTEGER
  typedef int            INTEGER;
  #define INTEGER_DEFINED
#endif
#ifndef UINTEGER_DEFINED
  #undef UINTEGER
  typedef unsigned int   UINTEGER;
  #define UINTEGER_DEFINED
#endif
#ifndef ULINTEGER_DEFINED
  #undef ULINTEGER
  typedef unsigned long  ULINTEGER;
  #define ULINTEGER_DEFINED
#endif
#ifdef __MACH__
  #ifndef UINT64_DEFINED
    #undef UINT64
    typedef uint64_t     UINT64;
    #define UINT64_DEFINED
  #endif
#endif
#ifndef REAL_DEFINED
  #undef REAL
  typedef double         REAL;
  #define REAL_DEFINED
#endif
#ifndef REAL_EXT_DEFINED
  #undef REAL_EXT
  typedef long double    REAL_EXT;
  #define REAL_EXT_DEFINED
#endif

#ifndef STRING_DEFINED
  #if defined(__linux__) || defined(__MACH__)
    typedef std::string    STRING;
  #elif defined(_WIN32)
    #ifdef _UNICODE
      typedef std::wstring STRING;
    #else
      typedef std::string  STRING;
    #endif
  #endif
  #define STRING_DEFINED
#endif

typedef std::basic_stringstream< STRING::value_type > STRINGSTREAM;
typedef std::basic_ostream< STRING::value_type >      OSTREAM;
typedef std::basic_istream< STRING::value_type >      ISTREAM;

typedef std::basic_streambuf< STRING::value_type >    STREAMBUFFER;
typedef std::basic_filebuf< STRING::value_type >      FILESTREAMBUFFER;

typedef boost::basic_format< STRING::value_type >     BOOSTFORMAT;

// getpid
#if defined(_WIN32)
  #define GETPID _getpid
#else // Linux, Mac & others
  #define GETPID getpid
#endif

///////////////////////////////////////////////////////////////
// Types specific to pSSAlib
namespace pssalib
{
  namespace datamodel
  {
    class DataModel;
  }

  /**
   * @param done Current sample index
   * @param total Total number of samples
   * @param percent Progress on current sample, in percent
   * @param ptrUser Pointer to a data structure supplied by user
   */
  typedef void (*FCN_REPORTPROGRESS_CALLBACK) (UINTEGER done, UINTEGER total, SHORT percent, void * ptrUser);

  /**
   * @param ptrDM Pointer to a \c DataModel object
   * @param time Current simulation time-point
   * @param ptrUser Pointer to a data structure supplied by user
   */
  typedef void (*FCN_REACTION_CALLBACK) (pssalib::datamodel::DataModel * ptrDM, REAL time, void * ptrUser);

  /**
   * @param ptrDM Pointer to a \c DataModel object
   * @param arPop Array to be filled with initial species' amounts
   * @param ptrUser Pointer to a data structure supplied by user
   */
  typedef void (*FCN_POPULATION_INITIALIZER) (pssalib::datamodel::DataModel * ptrDM, UINTEGER ** arPop, void * ptrUser);
}

///////////////////////////////////////////////////////////////
// Global functions specific to pSSAlib
namespace pssalib
{
  // None
} /* close namespace pssalib */

#endif /* PSSALIB_TYPEDEFS_H_ */
