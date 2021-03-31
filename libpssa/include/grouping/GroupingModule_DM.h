/**
 * @file GroupingModule_DM.h
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
 * Grouping module definition for the Gillespie's Direct Method
 */

#ifndef PSSALIB_GROUPING_GROUPINGMODULE_DM_H_
#define PSSALIB_GROUPING_GROUPINGMODULE_DM_H_

#include "./GroupingModule.h"

namespace pssalib
{
namespace grouping
{
  /**
   * \class GroupingModule_DM
   * \brief Fill in the datastructures for the Gillespie's Direct Method.
   * 
   * \copydetails GroupingModule
   */
  class GroupingModule_DM : public GroupingModule
  {
  ////////////////////////////////
  // Constructors
  public:
    // Default constructor
    GroupingModule_DM();

    // Copy constructor
    GroupingModule_DM(GroupingModule &);

    // Destructor
    virtual ~GroupingModule_DM();

  ////////////////////////////////
  // Methods
  public:
    // Calculate total propensity
virtual bool initialize(pssalib::datamodel::SimulationInfo *);
  };

}  } // close namespaces pssalib and grouping

#endif /* PSSALIB_GROUPING_GROUPINGMODULE_DM_H_ */
