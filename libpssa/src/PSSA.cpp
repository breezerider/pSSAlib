/**
 * @file PSSA.cpp
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
 * Implementation of the library's public interface
 */

#include "../include/PSSA.h"
#include "../include/datamodel/DataModel.h"
#include "../include/datamodel/DataModel_DM.h"
#include "../include/datamodel/DataModel_PDM.h"
#include "../include/datamodel/DataModel_SPDM.h"
#include "../include/datamodel/DataModel_PSSACR.h"

#include "../include/grouping/GroupingModule.h"
#include "../include/grouping/GroupingModule_DM.h"
#include "../include/grouping/GroupingModule_PDM.h"
#include "../include/grouping/GroupingModule_SPDM.h"
#include "../include/grouping/GroupingModule_PSSACR.h"

#include "../include/sampling/SamplingModule.h"
#include "../include/sampling/SamplingModule_DM.h"
#include "../include/sampling/SamplingModule_PDM.h"
#include "../include/sampling/SamplingModule_SPDM.h"
#include "../include/sampling/SamplingModule_PSSACR.h"

#include "../include/update/UpdateModule.h"
#include "../include/update/UpdateModule_DM.h"
#include "../include/update/UpdateModule_PDM.h"
#include "../include/update/UpdateModule_SPDM.h"
#include "../include/update/UpdateModule_PSSACR.h"

#include "../include/datamodel/SimulationInfo.h"

#include "../include/util/FileSystem.h"
#include "../include/util/Timing.h"

// output for all of them
#ifdef PSSA_MODULE_LABEL
#  undef  PSSA_MODULE_LABEL
#  define PSSA_MODULE_LABEL pssalib::datamodel::SimulationInfo::ofNone
#endif

namespace pssalib
{
  typedef struct tagTimingInfo
  {
    REAL t;     //!< Time spent simulating the trial
    UINTEGER n; //!< Number of reactions fired
  } TimingInfo;

  ///////////////////////////////
  // Constructors

  // Default constructor
  PSSA::PSSA()
    : ptrProgrCallback(NULL)
    , ptrReactionCallback(NULL)
    , ptrProgrCallbackUserData(NULL)
    , ptrReactionCallbackUserData(NULL)
    , ptrData(NULL)
    , ptrGrouping(NULL)
    , ptrSampling(NULL)
    , ptrUpdate(NULL)
    , m_Method(M_Invalid)
  {
    // Do nothing
  }

  // Destructor
  PSSA::~PSSA()
  {
    setMethod(M_Invalid);
  }

  ///////////////////////////////
  // Getter & Setters

  /**
   * Converts an SSA id into a human-redable string.
   * @param m Simulation algorithm id
   * @return A string representation of a given SSA.
   */
  STRING PSSA::getMethodName(const PSSA::EMethod m)
  {
    switch(m) {
      case M_Invalid:return STRING("Invalid value");
      case M_DM:     return STRING("DM");
      case M_PDM:    return STRING("PDM");
      case M_PSSACR: return STRING("PSSACR");
      case M_SPDM:   return STRING("SPDM");
      default:       return STRING("Unknown method");
    }
    return STRING();
  }

  /**
   * Return an SSA id matching a human-readable name.
   * @param s Simulation algorithm name
   * @return A numeric value associated with a given SSA.
   */
  PSSA::EMethod PSSA::getMethodID(const STRING &s)
  {
    if((0 == s.compare(0,2,"dm"))||
       (0 == s.compare(0,13,"direct method"))||
       (0 == s.compare(0,27,"gillespie's direct method")))
    {
      return M_DM;
    }
    else if((0 == s.compare(0,3,"pdm"))||
            (0 == s.compare(0,32,"partial-propensity direct method")))
    {
      return M_PDM;
    }
    else if((0 == s.compare(0,6,"pssacr"))||
            (0 == s.compare(0,58,"partial-propensity ssa with composition-rejection sampling"))) // Partial-Propensity SSA with Composition-Rejection sampling
    {
      return M_PSSACR;
    }
    else if((0 == s.compare(0,4,"spdm"))||
            (0 == s.compare(0,40,"sorting partial-propensity direct method")))
    {
      return M_SPDM;
    }
    else
    {
      return M_Invalid;
    }
  }

  /**
   * Sets the reaction callback
   * @param fcnProgress Pointer to a callback function that conforms with the 
   *                    FCN_REPORTPROGRESS_CALLBACK prototype
   * @param user        Pointer to user data passed on as a callback argument
   */
  void PSSA::SetProgressCallback(FCN_REPORTPROGRESS_CALLBACK fcnProgress, void* user)
  {
    ptrProgrCallback = fcnProgress;
    ptrProgrCallbackUserData = user;
  }

  /**
   * Sets the reaction callback
   * @param fcnReaction Pointer to a callback function that conforms with the 
   *                    FCN_REACTION_CALLBACK prototype
   * @param user        Pointer to user data passed on as a callback argument
   */
  void PSSA::SetReactionCallback(FCN_REACTION_CALLBACK fcnReaction, void* user)
  {
    ptrReactionCallback = fcnReaction;
    ptrReactionCallbackUserData = user;
  }
 
  /**
   * Sets the simulation method
   * @param NewMethod New method to use, must be one of the \b Method enumeration values
   * @return \b true if the simulation engine was set up properly, \b false otherwise
   */
  bool PSSA::setMethod(EMethod NewMethod)
  {
#ifdef BOOST_NO_CXX11_SMART_PTR
    std::auto_ptr<datamodel::DataModel> tempData, oldData;
    std::auto_ptr<grouping::GroupingModule> tempGrouping;
    std::auto_ptr<sampling::SamplingModule> tempSampling;
    std::auto_ptr<update::UpdateModule> tempUpdate;
#else
    std::unique_ptr<datamodel::DataModel> tempData, oldData;
    std::unique_ptr<grouping::GroupingModule> tempGrouping;
    std::unique_ptr<sampling::SamplingModule> tempSampling;
    std::unique_ptr<update::UpdateModule> tempUpdate;
#endif

    if(m_Method != NewMethod)
    {
      try
      {
        // Invalidate the object
        m_Method = M_Invalid;

        oldData.reset(ptrData);
        ptrData = NULL;

        delete ptrGrouping;
        ptrGrouping = NULL;

        delete ptrSampling;
        ptrSampling = NULL;

        delete ptrUpdate;
        ptrUpdate = NULL;

        switch(NewMethod)
        {
        // (Delayed) Gillespie's DM
        case M_DM:
          tempData.reset(new datamodel::DataModel_DM());
          tempGrouping.reset(new grouping::GroupingModule_DM());
          tempSampling.reset(new sampling::SamplingModule_DM());
          tempUpdate.reset(new update::UpdateModule_DM());
          break;
        // (Delayed) Partial Propensity Direct Method
        case M_PDM:
          tempData.reset(new datamodel::DataModel_PDM());
          tempGrouping.reset(new grouping::GroupingModule_PDM());
          tempSampling.reset(new sampling::SamplingModule_PDM());
          tempUpdate.reset(new update::UpdateModule_PDM());
          break;
        // (Delayed) PSSA with Composition-Rejection Sampling
        case M_PSSACR:
          tempData.reset(new datamodel::DataModel_PSSACR());
          tempGrouping.reset(new grouping::GroupingModule_PSSACR());
          tempSampling.reset(new sampling::SamplingModule_PSSACR());
          tempUpdate.reset(new update::UpdateModule_PSSACR());
          break;
        // (Delayed) Sorting Partial Propensity Direct Method
        case M_SPDM:
          tempData.reset(new datamodel::DataModel_SPDM());
          tempGrouping.reset(new grouping::GroupingModule_SPDM());
          tempSampling.reset(new sampling::SamplingModule_SPDM());
          tempUpdate.reset(new update::UpdateModule_SPDM());
          break;
        // Unset
        case M_Invalid:
        // Illegal parameter value
        default:
          NewMethod = M_Invalid;
        }

        // swap the data strucutres
        if((M_Invalid != NewMethod)&&
           (NULL != oldData.get())&&
           (NULL != tempData.get()))
          tempData.get()->swap(*oldData.get());

        // Store new value
        m_Method = NewMethod;
      }
      catch (std::bad_alloc & e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: " << e.what() << std::endl;
        m_Method = M_Invalid;
      }
      catch (std::runtime_error & e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error: " << e.what() << std::endl;
        m_Method = M_Invalid;
      }

      if(M_Invalid != m_Method)
      {
        ptrData = tempData.release();
        ptrGrouping = tempGrouping.release();
        ptrSampling = tempSampling.release();
        ptrUpdate = tempUpdate.release();
        return true;
      }
      else
        return false;
    }
    else
      return (M_Invalid != m_Method);
  }

  ///////////////////////////////////////////////////
  // Simulation

  bool PSSA::initSimulation(datamodel::SimulationInfo* ptrSimInfo)
  {
    //////////////////////////////
    // Attach SimulationInfo object to this instance
    if(!ptrSimInfo->attachPSSA(this))
    {
      // Output setup failed, do not attempt to write out
      return false;
    }

    //////////////////////////////
    // Self-check
    if(!isValid())
    {
      PSSA_ERROR(ptrSimInfo, << "the pSSAlib engine has not been initialized properly.\n");
      return false;
    }

    // Check input params
    if(!ptrSimInfo->isValid())
    {
      PSSA_ERROR(ptrSimInfo, << "simulation parameters are invalid.\n");
      return false;
    }

    // Reset timers
    ptrSimInfo->dTimeCheckpoint = 0.0;
    ptrSimInfo->dTimeSimulation = 0.0;

    return true;
  }

  bool PSSA::deinitSimulation(datamodel::SimulationInfo* ptrSimInfo)
  {
    // Dettach SimulationInfo object from this instance
    ptrSimInfo->detachPSSA();

    return true;
  }

  ///////////////////////////////
  // Protected members

  bool PSSA::setupForSampling(datamodel::SimulationInfo* ptrSimInfo)
  {
    //////////////////////////////
    // Initialize the simulation engine
    if(!initSimulation(ptrSimInfo))
    {
      PSSA_ERROR(ptrSimInfo, << "failed to initialize the simulation engine for sampling.\n");
      return false;
    }

    //////////////////////////////
    // Initialize the data structures
    if(!ptrGrouping->preinitialize(ptrSimInfo))
    {
      PSSA_ERROR(ptrSimInfo, << "failed to initialize data structures.\n");
      return false;
    }

    //////////////////////////////
    // Process user settings
    if(!ptrSimInfo->processSettings())
    {
      PSSA_ERROR(ptrSimInfo, << "failed to process user settings.\n");
      return false;
    }

    return true;
  }

  bool PSSA::runSamplingLoop(datamodel::SimulationInfo* ptrSimInfo)
  {
    static std::stringstream ssTemp;

    boost::scoped_array<TimingInfo> arTiming(NULL);
    boost::scoped_array<UINTEGER> arFinalPops(NULL);
    UINTEGER *ptrarFinalPops;

    PSSA_INFO(ptrSimInfo, << "# of species ids in simulation output "
      << ptrSimInfo->m_arSpeciesIdx.size() << ".\n");

    // Species names
    if(ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofSpeciesIDs))
    {
#ifdef HAVE_MPI
      bool bSpeciesNamesOK = true;
#endif
      if((ptrSimInfo->m_arSpeciesIdx.size() > 0)
#ifdef HAVE_MPI
        &&getMPIWrapperInstance().isMaster()
#endif
        )
      {
        PSSA_INFO(ptrSimInfo, << "writing species ids to a stream.\n");

        std::ostream & osLocal = ptrSimInfo->getOutputStream(datamodel::SimulationInfo::ofSpeciesIDs);
        if(!osLocal.good())
        {
          PSSA_ERROR(ptrSimInfo, << "invalid stream state!\n");
          ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofSpeciesIDs);
#ifndef HAVE_MPI
          return false;
#else
          bSpeciesNamesOK = false;
#endif
        }
        else
        {
          for(UINTEGER i = 0; i < ptrSimInfo->m_arSpeciesIdx.size(); ++i)
            osLocal << ptrData->getSpecies(ptrSimInfo->m_arSpeciesIdx[i])->getId() 
              << std::endl;

          ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofSpeciesIDs);

          PSSA_INFO(ptrSimInfo, << "species ids written to stream.\n");
        }
      }
#ifdef HAVE_MPI
      bSpeciesNamesOK = getMPIWrapperInstance().sync_results(bSpeciesNamesOK);
      if(!bSpeciesNamesOK)
        return false;
#endif
    }

    // Time points
#ifdef HAVE_MPI
    bool bTimePointsOK = true;
#endif
    if((ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofTimePoints))
#ifdef HAVE_MPI
       &&getMPIWrapperInstance().isMaster()
#endif
      )
    {
      PSSA_INFO(ptrSimInfo, << "writing time points to a stream.\n");

      std::ostream & osLocal = ptrSimInfo->getOutputStream(datamodel::SimulationInfo::ofTimePoints);
      if(!osLocal.good())
      {
          PSSA_ERROR(ptrSimInfo, << "invalid stream state!\n");
          ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofTimePoints);
#ifndef HAVE_MPI
        return false;
#else
        bTimePointsOK = false;
#endif
      }
      else
      {
        // Output to file
        UINTEGER unNumTimePoints = timing::getNumTimePoints(ptrSimInfo->dTimeStart,
                                                            ptrSimInfo->dTimeEnd,
                                                            ptrSimInfo->dTimeStep);
        REAL dTP = ptrSimInfo->dTimeStart;
        for(UINTEGER n = 0; n < (unNumTimePoints - 1); n++, dTP += ptrSimInfo->dTimeStep)
        {
          osLocal << dTP << "\n";
        }
        osLocal << ptrSimInfo->dTimeEnd << "\n";

        ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofTimePoints);

        PSSA_INFO(ptrSimInfo, << "time points written to stream.\n");
      }
    }
#ifdef HAVE_MPI
    bTimePointsOK = getMPIWrapperInstance().sync_results(bTimePointsOK);
    if(!bTimePointsOK)
      return false;
#endif

    // Initialize the par-for loop
#ifdef HAVE_MPI
    getMPIWrapperInstance().pre_spread(ptrSimInfo);//, unSamples);
#endif

    bool bAllOK = true;
    try
    {
      // Timing
      if (ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofTiming))
      {
        arTiming.reset(
#ifdef HAVE_MPI
          (TimingInfo*)getMPIWrapperInstance().spread_alloc(ptrSimInfo, sizeof(TimingInfo))
#else
          new TimingInfo[ptrSimInfo->unSamplesTotal]
#endif
                );
      }

      // Species populations at the final time point
      if (ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofFinalPops))
      {
        arFinalPops.reset(
#ifdef HAVE_MPI
          (UINTEGER*)getMPIWrapperInstance().
            spread_alloc(ptrSimInfo,
                        sizeof(UINTEGER)*ptrData->getSubvolumesCount()*ptrSimInfo->m_arSpeciesIdx.size())
#else
          new UINTEGER[ptrSimInfo->unSamplesTotal*ptrData->getSubvolumesCount()*ptrSimInfo->m_arSpeciesIdx.size()]
#endif
                );
        ptrarFinalPops = arFinalPops.get();
      } else if (ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofRawFinalPops)) {
        ptrarFinalPops = ptrSimInfo->ptrarRawPopulations;
      }
    }
    catch(std::bad_alloc & e)
    {
      PSSA_ERROR(ptrSimInfo, << e.what() << ": unable to allocate memory.\n");
      bAllOK = false;
    }
#ifdef HAVE_MPI
    bAllOK = getMPIWrapperInstance().sync_results(bAllOK);
#endif
    if(!bAllOK)
      return false;


    UINTEGER n = 0,  n_it = 0;
#ifdef HAVE_MPI
#error "disable building the version with MPI support. Just remove this to reenable"
    while(getMPIWrapperInstance().spread(ptrSimInfo,n))
#else
    for(; n < ptrSimInfo->unSamplesTotal; ++n)
#endif
    {
      // Clear streams
      ssTemp.str(STRING());

      // Initialize data structures
      bool bResult = ptrGrouping->initialize(ptrSimInfo);
#ifdef HAVE_MPI
      bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
      if(!bResult)
      {
        // Failed, report & exit
        PSSA_ERROR(ptrSimInfo, << "failed to initialize data structures.\n");
        return false;
      }
      // Optional post-initialisation step
      ptrGrouping->postInitialize(ptrSimInfo);

      // Timing
      REAL tTrial = 0.0;
      UINTEGER unReactions = 0;

      // Start timing
      bResult = ptrSimInfo->beginTrial(n);
#ifdef HAVE_MPI
      bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
      if(!bResult)
      {
        // Failed, report & exit
        PSSA_ERROR(ptrSimInfo, << "failed to initialize timing.\n");
        return false;
      }

      /////////////////////////////////
      // Run the internal loop
      while(ptrSimInfo->isRunning())
      {
        bool bSimResult = ptrSampling->getSample(ptrSimInfo);

        ptrSimInfo->doOutput();

        if(bSimResult)
        {
          bSimResult = ptrUpdate->doUpdate(ptrSimInfo);
          ++unReactions;
        }
        else
        {
            PSSA_WARNING(ptrSimInfo, << "sampling step failed after "
              << unReactions << " reactions.\n"
              << "sample = " << n << "\nsimulation time = "
              << ptrSimInfo->dTimeSimulation << "\ntotal propensity = "
              << ptrData->dTotalPropensity << std::endl);
            // is it an absorbing state?
            if(isinf(ptrSimInfo->dTimeSimulation))
              break; // exit the loop silently
            else
            {
              bResult = false; // fail
              break;
            }
        }

        if(bSimResult)
        {
          if (NULL != ptrReactionCallback)
          {
            ptrReactionCallback(ptrData,
              ptrSimInfo->dTimeSimulation,
              ptrReactionCallbackUserData);
          }
        }
        else // update failed
        {
          // is the simulation still running?
          if(!ptrSimInfo->isRunning())
            break; // exit the loop silently
          else
          {
            bResult = false; // fail
            break;
          }
        }
        PSSA_TRACE(ptrSimInfo, << "reaction " << unReactions
          << " simulation time = "
          << ptrSimInfo->dTimeSimulation << "; total propensity = "
          << ptrData->dTotalPropensity << std::endl);

        // handle external interruption
        if(ptrSimInfo->bInterruptRequested)
        {
          bResult = false;
          break;
        }
      }

      // End timing
      if(bResult)
        tTrial = ptrSimInfo->endTrial();
      else
        tTrial = 0.0;

#ifdef HAVE_MPI
      bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
      if(!bResult)
      {
        PSSA_ERROR(ptrSimInfo, << "simulation terminated unexpectedly!\n" 
          << "Sample : " << n << "\nSimulation time :"
          << ptrSimInfo->dTimeSimulation << "\nTotal propensity: "
          << ptrData->dTotalPropensity << "\nPrevious reaction : "
          << ptrData->getReactionWrapper(ptrData->mu).toString() << std::endl);
        return false; // fail
      }

      ////////////////////////
      // Store simulation results

      // Store the population at final time point
      if(ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofFinalPops)||ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofRawFinalPops))
      {
        for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); svi++) {
          datamodel::detail::Subvolume & subvol = ptrData->getSubvolume(svi);
          for(UINTEGER i = 0; i < ptrSimInfo->m_arSpeciesIdx.size(); i++) {
            ptrarFinalPops[(n_it*ptrData->getSubvolumesCount() + svi)*ptrSimInfo->m_arSpeciesIdx.size() + i] =
              subvol.population(ptrSimInfo->m_arSpeciesIdx[i]);
          }
        }
      }
      else
        PSSA_INFO(ptrSimInfo, << "final populations are not collected.\n");

      // Store the timing information
      if (ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofTiming))
      {
        PSSA_INFO(ptrSimInfo, << "timing info at iteration " << n_it
                  << ": \ttime = " << tTrial << "; NumReactions = "
                  << unReactions << std::endl);
        arTiming[n_it].t = tTrial;
        arTiming[n_it].n = unReactions;
      }
      else
        PSSA_INFO(ptrSimInfo, << "timing information is not collected.\n");

      n_it++;
    }

    bool bResult = true;
#ifdef HAVE_MPI
    if(getMPIWrapperInstance().post_spread(ptrSimInfo))
    {
      PSSA_TRACE(ptrSimInfo, << "Waiting for other processes to finish their task...\n");

      // Initialize data structures
      bResult = getMPIWrapperInstance().sync_results(bResult);
      if(!bResult)
        return false;

      PSSA_TRACE(ptrSimInfo, << "sync after GroupingModule::initialize()\n");

      // Start timing
      bResult = getMPIWrapperInstance().sync_results(bResult);
      if(!bResult)
        return false;

      PSSA_TRACE(ptrSimInfo, << "sync before the simulation loop\n");

      // After the simulation loop
      bResult = getMPIWrapperInstance().sync_results(bResult);
      if(!bResult)
        return false;

      PSSA_TRACE(ptrSimInfo, << "sync after the simulation loop\n");

      PSSA_TRACE(ptrSimInfo, << "Other processes' tasks are finished!\n");
    }
    else
    {
      PSSA_TRACE(ptrSimInfo, << "post_spread returned false.\n");
    }
#endif

    // Timing
    if(ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofTiming))
    {
      PSSA_INFO(ptrSimInfo, << "collecting timing info\n");
      TimingInfo * arCumTiming = NULL;
#ifdef HAVE_MPI
      bool bTimingOK = true;
      if(getMPIWrapperInstance().spread_collect(ptrSimInfo, arTiming.get(), (void **)&arCumTiming, sizeof(TimingInfo)))
      {
        if(getMPIWrapperInstance().isMaster())
        {
#else
          arCumTiming = arTiming.get();
#endif
          std::ostream & osLocal = ptrSimInfo->getOutputStream(datamodel::SimulationInfo::ofTiming);
          if(!osLocal.good())
          {
            PSSA_ERROR(ptrSimInfo, << "timing stream is invalid!\n");
#ifndef HAVE_MPI
            return false;
#else
            bTimingOK = false;
#endif
          }
          else
          {
            for(UINTEGER i = 0; i < ptrSimInfo->unSamplesTotal; i++)
              osLocal << arCumTiming[i].t << PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER
                << arCumTiming[i].n << std::endl;

            ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofTiming);
            PSSA_INFO(ptrSimInfo, << "timing info written to stream.\n");
          }
#ifdef HAVE_MPI
          // clean-up
          delete [] arCumTiming;
        }
      }
      else
      {
        bTimingOK = false;
      }
      bTimingOK = getMPIWrapperInstance().sync_results(bTimingOK);
      if(!bTimingOK)
        return false;
#endif
    }

    // Final time point populations
    if(ptrSimInfo->isLoggingOn(datamodel::SimulationInfo::ofFinalPops))
    {
      PSSA_INFO(ptrSimInfo, << "collecting final populations\n");
      UINTEGER *arCumFinalPops = NULL;
#ifdef HAVE_MPI
      bool bFinalPopsOK = true;
      if(getMPIWrapperInstance().spread_collect(ptrSimInfo, ptrarFinalPops,
        (void **)&arCumFinalPops, sizeof(UINTEGER)*ptrData->getSubvolumesCount()*ptrSimInfo->m_arSpeciesIdx.size()))
      {
        if(getMPIWrapperInstance().isMaster())
        {
#else
          arCumFinalPops = ptrarFinalPops;
#endif
          std::ostream & osLocal = ptrSimInfo->getOutputStream(datamodel::SimulationInfo::ofFinalPops);
          if(!osLocal.good())
          {
            PSSA_ERROR(ptrSimInfo, << "final populations stream is invalid!\n");
#ifndef HAVE_MPI
            return false;
#else
            bFinalPopsOK = false;
#endif
          }
          else
          {
            for(UINTEGER k = 0; k < ptrSimInfo->unSamplesTotal; k++)
            {
              for(UINTEGER si = 0; si < ptrData->getSubvolumesCount(); si++)
              {
                if (si != 0) osLocal << PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER;
                for(UINTEGER i = 0; i < ptrSimInfo->m_arSpeciesIdx.size(); i++)
                {
                  osLocal << arCumFinalPops[si*ptrSimInfo->m_arSpeciesIdx.size()
                    + i + k*ptrData->getSubvolumesCount()*ptrSimInfo->m_arSpeciesIdx.size()]
                    << PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER;
                }
              }
              osLocal << '\n';
            }
            ptrSimInfo->resetOutputStream(datamodel::SimulationInfo::ofFinalPops);
            PSSA_INFO(ptrSimInfo, << "final populations written to stream.\n");
          }
#ifdef HAVE_MPI
          // clean-up
          delete [] arCumFinalPops;
        }
      }
      else
      {
        bFinalPopsOK = false;
      }

      bFinalPopsOK = getMPIWrapperInstance().sync_results(bFinalPopsOK);
      if(!bFinalPopsOK)
        return false;
#endif
    }

    PSSA_INFO(ptrSimInfo, << "Sampling successfully completed, total iterations "
      << n_it << "; last sample #" << n << "; total samples "
      << ptrSimInfo->unSamplesTotal << std::endl);

    return bResult;
  }

  /**
   * This function samples \c ptrSimInfo->unSamplesTotal trajectories and outputs them to a series of files, 
   * starting at \a time \c = \c ptrSimInfo->dTimeStart to \a time \c = \c ptrSimInfo->dTimeEnd seconds and saving 
   * the output to \c ptrSimInfo->strOutput every \c ptrSimInfo->dTimeStep seconds.
   *
   * @param ptrSimInfo datamodel::SimulationInfo* Simulation information object associated with this run.
   *
   * @return logical \a "true" if all trials finished successfully,
   * \a "false" otherwise.
   *
   */
  bool PSSA::run(datamodel::SimulationInfo *ptrSimInfo)
  {
    bool bResult = setupForSampling(ptrSimInfo);
#ifdef HAVE_MPI
    bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
    if(!bResult)
      return false;

    bResult = runSamplingLoop(ptrSimInfo);

    // Clean-up
    deinitSimulation(ptrSimInfo);

    return bResult;
  }

  /**
   * This function samples given number of trajectories (\c ptrSimInfo->arSamples[0]) and outputs them to a series of files, 
   * starting at \a time \c = \c ptrSimInfo->dTimeStart to \a time \c = \c ptrSimInfo->dTimeEnd seconds and saving 
   * the output to \c ptrSimInfo->strOutput every \c ptrSimInfo->dTimeStep seconds.
   *
   * @param ptrSimInfo datamodel::SimulationInfo* Simulation information object associated with this run.
   *
   * @return logical \a "true" if all trials finished successfully,
   * \a "false" otherwise.
   *
   */
  bool PSSA::run_avg(datamodel::SimulationInfo *ptrSimInfo)
  {
    bool bResult = setupForSampling(ptrSimInfo);
#ifdef HAVE_MPI
    bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
    if(!bResult)
      return false;

    //////////////////////////////
    // Run the simulation
    UINTEGER prevOutputFalgs = ptrSimInfo->unOutputFlags;

    ptrSimInfo->unOutputFlags |= 
      datamodel::SimulationInfo::ofTrajectory |
      datamodel::SimulationInfo::ofTimePoints;
    ptrSimInfo->unOutputFlags &= ~(
      datamodel::SimulationInfo::ofFinalPops |
      datamodel::SimulationInfo::ofTiming
    );

    bResult = runSamplingLoop(ptrSimInfo);

    //
    // Clean-up
    deinitSimulation(ptrSimInfo);

    // restore previous output flags
    ptrSimInfo->unOutputFlags = prevOutputFalgs;

    return bResult;
  }

  /**
   * This function simulates given in \c simInfo->arSamples number of trajectories and outputs the population at \a time \c = \c simInfo->dTimeEnd 
   * to a file named \c "PSSA_final_values.txt". For each element in \c arTrials a separate output is produced and the data from previous trials (if any) is
   * reused. At the same time it computes the amount of time spent one each run and the number of reactions fired and stores them to \c "PSSA_probability_timing.txt". 
   * Depending on the value of \c simInfo->stStatistics, it can produce the one-dimensional histograms for each of the species (\a ST_SINGLE) or cumulative histogram 
   * of all of the species (\a ST_MULTI). Only species in \c simInfo->pArSpeciesIds are considered.
   *
   * @param simInfo datamodel::SimulationInfo* Simulation information object associated with this run.
   *
   * @return logical @a "true" if the all trials finish successfully,
   * @a "false" otherwise.
   *
   */
  bool PSSA::run_hist(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    bool bResult = setupForSampling(ptrSimInfo);
#ifdef HAVE_MPI
    bResult = getMPIWrapperInstance().sync_results(bResult);
#endif
    if(!bResult)
      return false;

    //////////////////////////////
    // Run the simulation
    UINTEGER prevOutputFalgs = ptrSimInfo->unOutputFlags;

    ptrSimInfo->unOutputFlags |= 
      datamodel::SimulationInfo::ofFinalPops |
      datamodel::SimulationInfo::ofTimePoints;
    ptrSimInfo->unOutputFlags &= ~(
      datamodel::SimulationInfo::ofTrajectory |
      datamodel::SimulationInfo::ofTiming
    );

    bResult = runSamplingLoop(ptrSimInfo);

    //
    // Clean-up
    deinitSimulation(ptrSimInfo);

    // restore previous output flags
    ptrSimInfo->unOutputFlags = prevOutputFalgs;

    return bResult;
  }
} // close namespace pssalib
