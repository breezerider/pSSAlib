/**
 * @file GroupingModule_SPDM.h
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
 * Grouping module definition for the Sorting Partial Propensity Direct Method
 */

#ifndef PSSALIB_GROUPING_GROUPINGMODULE_SPDM_H_
#define PSSALIB_GROUPING_GROUPINGMODULE_SPDM_H_

#include "./GroupingModule_PDM.h"

namespace pssalib
{
namespace grouping
{
  /**
   * @class GroupingModule_SPDM
   * @brief Fill in the datastructures for the Sorting Partial-Propensity 
   * Direct Method.
   * 
   * @copydetails GroupingModule
   */
  class GroupingModule_SPDM : public GroupingModule_PDM
  {
  ////////////////////////////////
  // Constructors
  public:
    // Default constructor
    GroupingModule_SPDM();

    // Copy constructor
    GroupingModule_SPDM(GroupingModule &);

    // Destructor
virtual ~GroupingModule_SPDM();

  ////////////////////////////////
  // Methods
  public:
    // Initialize data structures (called before each trial)
virtual bool initialize(pssalib::datamodel::SimulationInfo *);
  };

}  } // close namespaces pssalib and grouping

#endif /* PSSALIB_GROUPING_GROUPINGMODULE_SPDM_H_ */
