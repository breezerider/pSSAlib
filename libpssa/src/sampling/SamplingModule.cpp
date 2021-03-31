/**
 * @file SamplingModule.cpp
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
 * Implementation of the generic Sampling module:
 *   - sampling reaction time & number using the Gillespie's method
 *   - composition-rejection sampling of the reactor subvolume
 */

#include "../../include/stdheaders.h"
#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/sampling/SamplingModule.h"

namespace pssalib
{
namespace sampling
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  SamplingModule::SamplingModule() 
    : m_ptrRNG(NULL)
  {
    gsl_rng_env_setup();
    const gsl_rng_type * ptrRNGtype = gsl_rng_default;
    m_ptrRNG = gsl_rng_alloc(ptrRNGtype);

#ifndef PSSALIB_ENGINE_CHECK
    if(NULL == getenv("GSL_RNG_SEED"))
    {
      unsigned long seed = 1;
#ifdef HAVE_MPI
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      seed *= boost::hash_value<int>(++rank);
#endif
      seed *= boost::hash_value<std::clock_t>(std::clock());

      gsl_rng_set(m_ptrRNG, seed);
    }
#endif
  }

  //! Destructor
  SamplingModule::~SamplingModule()
  {
    if(m_ptrRNG)
    {
      gsl_rng_free(m_ptrRNG);
      m_ptrRNG = NULL;
    }
  }

  ////////////////////////////////
  // Methods

  //! Fill in the datastructure with random samples
  bool SamplingModule::getSample(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

    // Sample time
    if(!sampleTime(ptrSimInfo))
    {
      PSSA_ERROR(ptrSimInfo, << "could not sample next reaction time!\n");
      return false;
    }

    // Sample volume
    if(0 != ptrData->getDimsCount())
    {
      if(!sampleVolume(ptrSimInfo))
      {
        PSSA_ERROR(ptrSimInfo, << "could not sample next reaction volume!\n");
        return false;
      }
    }

    // Sample reaction
    if(!sampleReaction(ptrSimInfo))
    {
      PSSA_ERROR(ptrSimInfo, << "could not sample next reaction index!\n");
      return false;
    }
    else
    {
      PSSA_TRACE(ptrSimInfo, << "sampled reaction #" << ptrData->mu 
        << " : " << ptrData->getReactionWrapper(ptrData->mu).toString()
        << std::endl);
    }

    // Sample diffusion destination if necessary.
    if(ptrData->getReactionWrapper(ptrData->mu).isDiffusive())
    {
      // With isotropic diffusion we don't need a linear search.
      UINTEGER idxDestSubVol = (UINTEGER)
        (gsl_rng_uniform_pos(m_ptrRNG) * 2 * ptrData->getDimsCount());
      ptrData->nu_D = ptrData->getSubvolume(ptrData->nu).neighbour(idxDestSubVol);
      PSSA_TRACE(ptrSimInfo, << "sampled destination volume : source = " << ptrData->nu
        << "; destination = " << ptrData->nu_D << std::endl);
    }

    return true;
  }

  //! Sample the time interval till next reaction
  bool SamplingModule::sampleTime(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

    // Generate a random number
    REAL r;
    do
    {
      r = gsl_rng_uniform (m_ptrRNG);
    } while (r == 0);

    if(ptrData->vQueuedReactions.empty())
    {
      // check if we have reached an absorbing state
      if(ptrData->dTotalPropensity <= 0.0)
      {
        ptrSimInfo->dTimeSimulation = std::numeric_limits<REAL>::infinity();
        PSSA_WARNING(ptrSimInfo, << "zero or negative propensity ==> simulation reached an absorbing state.\n");
        return false; // We have reached an absorbing state - exit
      }
      else
        ptrSimInfo->dTimeSimulation -= log(r) / ptrData->dTotalPropensity;

      PSSA_TRACE(ptrSimInfo, << "sampled time = " << ptrSimInfo->dTimeSimulation << "; total propensity = " << ptrData->dTotalPropensity << std::endl);
    }
    else
    {
      REAL T1, T2, at, F;
      std::vector<pssalib::datamodel::DataModel::DelayedReaction>::iterator 
        curReaction = ptrData->vQueuedReactions.begin();

      // Initialize
      T1 = ptrSimInfo->dTimeSimulation;
      T2 = curReaction->time;
      at = ptrData->dTotalPropensity * (T2 - T1);
      F  = 1.0 - exp(-at);

      while(F < r)
      {
        // Reached the end of the simulation
        ptrData->mu = curReaction->index;
        ptrSimInfo->dTimeSimulation = T2;
        // Write to file & update
        if(!ptrSimInfo->UpdateCallback()) 
        {
          PSSA_WARNING(ptrSimInfo, << "sampling failed: could not perform a delayed update!\n");
          return false;
        }

        // Check if we're still within simulation timespan
        if(T2 > ptrSimInfo->dTimeEnd)
          return false;

        // Remove the delayed reaction we just fired
        ptrData->vQueuedReactions.erase(curReaction);

        T1 = T2;
        // If all delayed reaction have already been fired
        // just calculate the time-step and exit.
        if(ptrData->vQueuedReactions.empty())
        {
          ptrSimInfo->dTimeSimulation = T1 - (gsl_log1p(-r) + at) / ptrData->dTotalPropensity;
          return true;
        }

        curReaction =  ptrData->vQueuedReactions.begin();
        T2          =  curReaction->time;
        at          += ptrData->dTotalPropensity * (T2 - T1);
        F           =  1.0 - exp(-at);
      }
      ptrSimInfo->dTimeSimulation = T2 - (gsl_log1p(-r) + at) / ptrData->dTotalPropensity;

      PSSA_TRACE(ptrSimInfo, << "sampled time (with delays) = "
        << ptrSimInfo->dTimeSimulation << std::endl);
    }

    return true;
  }

  //! Sample next reactor index using pSSA-CR
  bool SamplingModule::sampleVolume(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel * ptrData = ptrSimInfo->getDataModel();

    bool success;
//     if (ptrData->unSubvolumes <= 1)
//     {
//       ptrData->nu = 0;
//       success = true;
//     } else {
      UINTEGER i;
      REAL r;
      success = crVolumeSampler.Sample(&ptrData->crsdVolume, m_ptrRNG, 
                                       ptrData->dTotalPropensity, i, r);
      PSSA_TRACE(ptrSimInfo, << "sampled reactor subvolume = [ i=" << i << "; r=" << r << "]\n");
      if (!success) {
        ptrData->nu = 0;
        {
          PSSA_ERROR(ptrSimInfo, << "sampling did not converge within "
            "pre-defined number of iterations.\n");
        }
      } else {
        ptrData->nu = i;
      }
//     }

    return success;
  }

  /**
   * Set the seed of the random number generator
   * @param seed New seed value
   * 
   */
  void SamplingModule::set_rng_seed(UINTEGER seed)
  {
    gsl_rng_set(m_ptrRNG, seed);
  }

}  } // close namespaces pssalib and sampling
