/**
 * @file PSSA.h
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
 * This file contains the definition of the API.
 * It provides main functions of the pSSA library:
 * - creation and initialisation of the simulation engine;
 * - running simulations using a selection of SSAs;
 * - two output modes for results
 *   - "run"-mode: outputs trajectories
 *   - "stat"-mode: outputs PDFs.
 */

#ifndef PSSALIB_PSSA_H_
#define PSSALIB_PSSA_H_

#include "typedefs.h"

#ifndef PSSALIB_FILENAME_LOG
#  ifdef HAVE_MPI
#  define PSSALIB_FILENAME_LOG "libpssa_%i.log"
#  else
#  define PSSALIB_FILENAME_LOG "libpssa.log"
#  endif
#endif

#ifndef PSSALIB_FILENAME_STATUS
#  ifdef HAVE_MPI
#  define PSSALIB_FILENAME_STATUS "status_%i.txt"
#  else
#  define PSSALIB_FILENAME_STATUS "status.txt"
#  endif
#endif

#ifndef PSSALIB_FILENAME_SPECIES_IDS
#define PSSALIB_FILENAME_SPECIES_IDS "species_ids.dat"
#endif

#ifndef PSSALIB_FILENAME_TIME_POINTS
#define PSSALIB_FILENAME_TIME_POINTS "time_points.dat"
#endif

#ifndef PSSALIB_FILENAME_TRAJECTORY
#define PSSALIB_FILENAME_TRAJECTORY "trajectory_%i.dat"
#endif

#ifndef PSSALIB_FILENAME_FINAL_TIME_POINT_POPULATIONS
#define PSSALIB_FILENAME_FINAL_TIME_POINT_POPULATIONS "populations.dat"
#endif

#ifndef PSSALIB_FILENAME_TIMING
#define PSSALIB_FILENAME_TIMING "timing.dat"
#endif

///////////////////////////////////
// Forward declarations
namespace pssalib
{
  namespace datamodel
  {
    class DataModel;
    class SimulationInfo;
  }

  namespace grouping
  {
    class GroupingModule;
  }

  namespace sampling
  {
    class SamplingModule;
  }

  namespace update
  {
    class UpdateModule;
  }
}

namespace pssalib
{
  ///////////////////////////////////
  //! Partial propensity stochastic reaction-diffusion simulation library
  /*!
   * Main simulation engine class.
   */
  class PSSA
  {
  /////////////////////////////////
  // Friend classes
  public:
    friend class datamodel::SimulationInfo;

  /////////////////////////////////
  // Data Types
  public:
    typedef enum tagEMethod
    {
      //! Invalid Method
      M_Invalid = 0x0000,
      //! Gillespie's Direct Method
      M_DM = 0x0001,
      //! Partial Propensity Direct Method
      M_PDM = 0x0002,
      //! pSSA with Composition-Rejection Sampling
      M_PSSACR = 0x0004,
      //! Sorting Partial Propensity Direct Method
      M_SPDM = 0x0008,
      //! All methods
      M_All  = 0x000F
    } EMethod;

  /////////////////////////////////
  // Attributes
  protected:
    //! Pointer to progress callback function
    FCN_REPORTPROGRESS_CALLBACK   ptrProgrCallback;
    //! Pointer to reaction callback function
    FCN_REACTION_CALLBACK         ptrReactionCallback;
    //! Pointer to progress callback user data
    void*                         ptrProgrCallbackUserData;
//     //! Pointer to error callback user data
//     void*                         ptrErrCallbackUserData;
    //! Pointer to reaction callback user data
    void*                         ptrReactionCallbackUserData;

    //! Pointer to a corresponding data model object
    datamodel::DataModel*         ptrData;
    //! Pointer to a corresponding grouping module object
    grouping::GroupingModule*     ptrGrouping;
    //! Pointer to a corresponding sampling module object
    sampling::SamplingModule*     ptrSampling;
    //! Pointer to a corresponding update module object
    update::UpdateModule*         ptrUpdate;

    //! ID of the simulation method
    EMethod                       m_Method;

  /////////////////////////////////
  // Constructors
  public:
    //! Default Constructor
    PSSA();

    //! Destructor
    ~PSSA();

    //! Run the simulation using custom settings
    bool run(datamodel::SimulationInfo * simInfo);
    //! Sample given number of trajectories, output them and compute statistics.
    bool run_avg(datamodel::SimulationInfo * simInfo);
    //! Compute histogram of species at a given time point over a given number of trials.
    bool run_hist(datamodel::SimulationInfo * simInfo);

  //////////////////////////////
  // Methods
  protected:

    //! Run self-checks & prepare the simInfo object for simulations
    bool initSimulation(datamodel::SimulationInfo* simInfo);

    //! Free resources allocated for the simulation
    bool deinitSimulation(datamodel::SimulationInfo* simInfo);

    //! Initialize the engine for sampling
    bool setupForSampling(datamodel::SimulationInfo* simInfo);

    //! Simulation driver
    bool runSamplingLoop(datamodel::SimulationInfo* simInfo);

  //////////////////////////////
  // Methods
  public:

    //! Returns current simulation SSA is
    inline EMethod getMethod() const { return m_Method; };

    //! Returns a human-readable name of an SSA id
    static STRING getMethodName(const EMethod m);

    //! Returns an SSA id matching the human-readable name provided
    static EMethod getMethodID(const STRING &s);

    //! Sets the method used in simulations
    bool setMethod(EMethod m);
    //! Sets the progress callback function and user data
    void SetProgressCallback(FCN_REPORTPROGRESS_CALLBACK fcnProgress, void* user);
    //! Sets the reaction callback function and user data
    void SetReactionCallback(FCN_REACTION_CALLBACK fcnReaction, void* user);

    /**
     * Getters for modules
     */

    //! Returns the name of DataModel object asociated with this instance of the engine
    inline STRING getModelName() const;
// 
//     //! Returns number of species in the model
//     inline UINTEGER getNumSpecies() const;

    //! Returns true if the engine was properly set up
    inline bool isValid() const
    {
      if(m_Method != M_Invalid)
        return (ptrData != NULL)&&
               (ptrGrouping != NULL)&&
               (ptrSampling != NULL)&&
               (ptrUpdate != NULL);
      return false;
    }
  };
}

// For the inline method (see below)
#include "./datamodel/DataModel.h"

namespace pssalib
{
  inline STRING PSSA::getModelName() const
  {
    return (isValid() ? ptrData->getName() : STRING());
  }
// 
//   inline UINTEGER PSSA::getNumSpecies() const
//   {
//     return ptrData->getSpeciesCount();
//   }
}

// To provide the SimulationInfo class
// just by including main header file
#include "./datamodel/SimulationInfo.h"

#endif /* PSSALIB_PSSA_H_ */
