/**
 * @file UpdateModule.cpp
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
 * Implementation of the generic Update module:
 *   - update simulation data structures with changes due to chemical reactions
 *   - update simulation data structures with changes due to diffusion reactions
 *   - schedule & perform of delayed updates
 */

#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/update/UpdateModule.h"

#include "../../include/util/Maths.h"

namespace pssalib
{
namespace update
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  UpdateModule::UpdateModule()
    : m_sriBegin(0)
    , m_sriEnd(0)
    , m_sriReactants(0)
  {
    // Do nothing
  }

  //! Destructor
  UpdateModule::~UpdateModule()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Schedule delayed reaction
  bool UpdateModule::scheduleDelayed(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

    // Store the reaction
    pssalib::datamodel::DataModel::DelayedReaction
      reaction(ptrData->mu, ptrSimInfo->dTimeSimulation +
               ptrData->getReactionWrapper(ptrData->mu).getDelay());
    ptrData->vQueuedReactions.push_back(reaction);
    std::sort(ptrData->vQueuedReactions.begin(), ptrData->vQueuedReactions.end());

    return true;
  }

  //! Determine update mode & perform update
  bool UpdateModule::doUpdate(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  { 
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel * ptrData =
      ptrSimInfo->getDataModel();
    m_ptrReactionWrapper = &(ptrData->getReactionWrapper(ptrData->mu));
    m_ptrSubvolumeSrc = &(ptrData->getSubvolume(ptrData->nu));

    REAL totalPropensityChange = m_ptrSubvolumeSrc->dTotalPropensity;
    bool bUpdateOK = true;

    // Update species population
    if(m_ptrReactionWrapper->isDiffusive())
    {
      m_ptrSubvolumeDst = &(ptrData->getSubvolume(ptrData->nu_D));

      totalPropensityChange += m_ptrSubvolumeDst->dTotalPropensity;

      // update population
      const UINTEGER index = m_ptrReactionWrapper->getSpecies()->getIndex();
      m_ptrSubvolumeSrc->population_update(index, -1);
      m_ptrSubvolumeDst->population_update(index,  1);

      // update method data structures
      bUpdateOK = updateSpeciesStructuresDiffusion(ptrSimInfo);

      totalPropensityChange -= m_ptrSubvolumeSrc->dTotalPropensity + m_ptrSubvolumeDst->dTotalPropensity;
    }
    else
    {
      m_sriReactants = m_ptrReactionWrapper->getReactantsCount();
      m_sriBegin = 0;
      m_sriEnd = m_ptrReactionWrapper->getSpeciesReferencesCount();

      if(m_ptrReactionWrapper->isSetDelay())
      {
        if(m_ptrReactionWrapper->isSetDelayConsuming())
        {
          if(ptrSimInfo->getDelayedUpdate())
            m_sriBegin = m_ptrReactionWrapper->getReactantsCount();
          else
            m_sriEnd = m_ptrReactionWrapper->getReactantsCount();
        }
        else if(!ptrSimInfo->getDelayedUpdate())
          return true; // nothing to update
      }

      // update population
      for(UINTEGER sri = m_sriBegin; sri < m_sriEnd; ++sri)
      {
        const pssalib::datamodel::detail::SpeciesReference * sr = 
          m_ptrReactionWrapper->getSpeciesReferenceAt(sri);

        if(sr->isConstant()) continue;
        m_ptrSubvolumeSrc->population_update(sr, (sri >= m_sriReactants));
      }

      // update method data structures
      bUpdateOK = updateSpeciesStructuresReaction(ptrSimInfo);

      totalPropensityChange -= m_ptrSubvolumeSrc->dTotalPropensity;
    }

    // Update global propensity
    ptrData->dTotalPropensity -= totalPropensityChange;

    if(!bUpdateOK)
    {
      PSSA_ERROR(ptrSimInfo, << "update failed: could not update subvolume structures." << std::endl);
      return false;
    }
    else
    {
      if(ptrSimInfo->isLoggingOn(pssalib::datamodel::SimulationInfo::ofTrace | PSSA_MODULE_LABEL))
      {
        STRINGSTREAM ssTemp;

        ssTemp << "update: " << (m_ptrReactionWrapper->isDiffusive() ? "pop src :" : "pop :");
        pssalib::datamodel::detail::Subvolume & subVol = ptrData->getSubvolume(ptrData->nu);
        for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
          ssTemp << " " << m_ptrSubvolumeSrc->population(si);
        ssTemp << "\t";

        // diffusion
        if(m_ptrReactionWrapper->isDiffusive())
        {

          ssTemp << "pop dest :";
          for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
            ssTemp << " " << m_ptrSubvolumeDst->population(si);
          ssTemp << "\t";
        }

        if(m_ptrReactionWrapper->isDiffusive())
          ssTemp << "propensity : src=" << m_ptrSubvolumeSrc->dTotalPropensity
            << "; dest=" << m_ptrSubvolumeDst->dTotalPropensity << ";  ";
        ssTemp << "tot prop=" << ptrData->dTotalPropensity;

        PSSA_TRACE(ptrSimInfo, << ssTemp.rdbuf() << std::endl);
      }
    }

    // Update compartment data structures
    if(ptrData->getSubvolumesCount() > 1)
    {
      if(!this->updateVolumeStructures(ptrSimInfo))
      {
        PSSA_ERROR(ptrSimInfo, << "update failed: could not update compartment structures." << std::endl);
        return false;
      }
    }

    return true;
  }

  //! Update volume data structures
  bool UpdateModule::updateVolumeStructures(pssalib::datamodel::SimulationInfo* ptrSimInfo)
  {
    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel * ptrData = ptrSimInfo->getDataModel();

    // Update source volume propensity.
    pssalib::datamodel::detail::Subvolume & subVol = ptrData->getSubvolume(ptrData->nu);

    UINTEGER k = (UINTEGER)std::floor(fabs(LOG2(subVol.dTotalPropensity / ptrData->crsdVolume.minValue))) + 1;

    ptrData->crsdVolume.updateValue(k, ptrData->nu, subVol.dTotalPropensity);

    if(ptrData->getReactionWrapper(ptrData->mu).isDiffusive())
    {
      // Update destination volume propensity.
      pssalib::datamodel::detail::Subvolume & subVol_D = ptrData->getSubvolume(ptrData->nu_D);
      //UINTEGER k = floor_log2((UINTEGER)(sv.dTotalPropensity / ptrData->crsdVolume.minValue)) + 1;
      k = (UINTEGER)std::floor(fabs(LOG2(subVol_D.dTotalPropensity / ptrData->crsdVolume.minValue))) + 1;
      ptrData->crsdVolume.updateValue(k, ptrData->nu_D, subVol_D.dTotalPropensity);
    }

    return true;
  }

}  } // close namespaces pssalib and update
