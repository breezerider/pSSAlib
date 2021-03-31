/**
 * @file GroupingModule.h
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
 * Grouping module parses the SBML model definition and translates
 * it into an internal method-specific definition in the data model
 */

#ifndef PSSALIB_GROUPING_GROUPINGMODULE_H_
#define PSSALIB_GROUPING_GROUPINGMODULE_H_

#include "../typedefs.h"
#include "../datamodel/DataModel.h"
#include "../datamodel/SimulationInfo.h"

#ifdef PSSA_MODULE_LABEL
  #undef PSSA_MODULE_LABEL
#endif
#define PSSA_MODULE_LABEL pssalib::datamodel::SimulationInfo::eofModuleGrouping

// Forward declarations
namespace pssalib
{
namespace datamodel
{
  class DataModel;
  class SimulationInfo;
} // close namespace datamodel

namespace grouping
{
  /**
   * @class GroupingModule
   * @brief This class serves as a base class for all grouping schemes 
   * offered by this library.
   * 
   * @details Translates the generic SBML model definition into 
   * a method-specific representation and fill in the relevant
   * data structures 
   * The respective methods of this class may be overloaded, offering a 
   * flexible way to extend its functionality.
   */
  class GroupingModule
  {
  ////////////////////////////////
  // Attributes
  protected:
    //! Flag for successful loading of data
    bool bDataLoaded;

  ////////////////////////////////
  // Constructors
  public:
    // Default constructor
    GroupingModule();

    // Copy constructor
    GroupingModule(GroupingModule &);

    // Destructor
    virtual ~GroupingModule();

  ////////////////////////////////
  // Methods
  public:
    // Parse the SBML model
  virtual bool preinitialize(pssalib::datamodel::SimulationInfo *);

    // Initialize simulation data structures (called before each trial)
  virtual bool initialize(pssalib::datamodel::SimulationInfo *);

    // Initialize composition-rejection sampler for subvolumes
  virtual void postInitialize(pssalib::datamodel::SimulationInfo *);
  };

}  } // close namespaces pssalib and grouping


#endif /* PSSALIB_GROUPING_GROUPINGMODULE_H_ */
