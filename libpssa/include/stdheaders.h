/**
 * @file stdheaders.h
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
 * Common header file with standard includes.
 */

#ifndef PSSALIB_STDHEADERS_H_
#define PSSALIB_STDHEADERS_H_

#ifdef HAVE_CONFIG_H
  #include "config.h" // current build configuration
#endif

#include <iostream> // Standard I/O streams
#include <fstream>  // Standard File I/O streams
#include <iomanip>  // Standard header for console output manipulation
#include <string.h> // Standard String header
#include <sstream>  // Standard String Stream Header
#include <stdarg.h> // Standard Variable Length Arguments Header

// time related headers
#if defined(__linux__)
  #include <time.h>
#elif defined(__MACH__)
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#elif defined(_WIN32)
  #define NOMINMAX
  #include <windows.h>
  #include <tchar.h>
#endif

// Standard C++ headers
#include <stdexcept>  // definition of std::exception
#include <cerrno>     // C-style errno macro

// definition of mkdir
#if defined(__linux__) || defined(__MACH__)
  #include <sys/stat.h>
  #include <sys/types.h>
#elif defined(_WIN32)
  #include <direct.h>
#endif

// libM
#include <limits>  // Need it for min & lowest

// GNU Scientific Library
#include <gsl/gsl_math.h>      // Maths constants
#include <gsl/gsl_rng.h>       // GSL random generators
#include <gsl/gsl_randist.h>   // Normally-distributed random numbers
#include <gsl/gsl_sf_gamma.h>  // Factorial
#include <gsl/gsl_const_num.h> // Avogadro's Number

// STL
#include <algorithm>      // STL algorithms
#include <vector>         // STL vector container
#include <memory>         // STL C-style pointer wrapper
#include <atomic>         // STL atomic definitions
#include <unordered_map>  // STL unordered map

// libSBML
#ifdef HAVE_LIBSBML
  #include <sbml/SBMLTypes.h> // Systems Biology Markup Language library
#endif

// MPI
#ifdef HAVE_MPI
  #include <mpi.h>
#endif

// Define hashmap type for PSSACR_Bins class
//#define __USE_GOOGLE_HASH_MAP

#ifdef __USE_GOOGLE_HASH_MAP
  #include <google/dense_hash_map>
#endif

// Boost headers
#include <boost/config.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
// #include <boost/unordered_map.hpp>
#include <boost/current_function.hpp>
// #include <boost/array.hpp>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/tee.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#endif /* PSSALIB_STDHEADERS_H_ */
