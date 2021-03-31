/**
 * @file Hreaction.hpp
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
 * Validation test case interface
 */

#include "typedefs.h"

#ifndef PSSALIB_VALIDATION_HREACTION_HPP_
#define PSSALIB_VALIDATION_HREACTION_HPP_

/**
 * @struct Hreaction
 * @brief Simple bimolecular reaction network
 */
struct Hreaction
{
  /**
   * Get test case name
   * 
   * @return A zero-terminated string representing the test case name
   */
  virtual STRING::const_pointer getName() const = 0;

  /**
   * Get ids of species to output
   * 
   * @return Vector of species ids for output
   */
  std::vector<STRING> * getSpeciesIds() const
  {
    std::vector<STRING> * ids = new std::vector<STRING>;
    ids->push_back(STRING("A"));

    return ids;
  }

  /**
   * Generate the SBML model
   * 
   * @return an @c SBMLDocument object representing the network
   */
  virtual LIBSBML_CPP_NAMESPACE::SBMLDocument * generateSBML() const = 0;

  /**
   * Compute the analytic PDF
   * 
   * @return an @c SBMLDocument object representing the network
   */
  virtual REAL computePDF(UINTEGER n) const = 0;
};

#endif /* PSSALIB_VALIDATION_HREACTION_HPP_ */
