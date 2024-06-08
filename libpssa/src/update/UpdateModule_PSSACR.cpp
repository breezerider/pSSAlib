/**
 * @file UpdateModule_PSSACR.cpp
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
 * Update module implementation for the Partial-Propensity SSA 
 * with Composition-Rejection sampling (Ramaswamy, 2010)
 */

#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/update/UpdateModule_PSSACR.h"

#include "../../include/util/Maths.h"

namespace pssalib
{
namespace update
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  UpdateModule_PSSACR::UpdateModule_PSSACR()
  {
    // Do nothing
  }

  //! Destructor
  UpdateModule_PSSACR::~UpdateModule_PSSACR()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods
  
  // Update per species data structures
  void UpdateModule_PSSACR::updateSpeciesStructures(pssalib::datamodel::SimulationInfo * ptrSimInfo,
                                                    pssalib::datamodel::DataModel_PSSACR * ptrPSSACRData,
                                                    pssalib::datamodel::detail::Subvolume_PSSACR & PSSACRSubVol,
                                                    const UINTEGER index)
  {
    bool bUpdateSelf = true;
    pssalib::datamodel::CompositionRejectionSamplerData & crsdSigma =
      PSSACRSubVol.crsdSigma;
    const REAL invMinSigma = 1.0 / crsdSigma.minValue;

    const UINTEGER adjusted_index = index + 1;

    PSSA_TRACE(ptrSimInfo,  << "updating reactions for species index " << adjusted_index << std::endl);

    for(UINTEGER l = 0, U3_rowlen = ptrPSSACRData->arU3.get_cols(adjusted_index); l < U3_rowlen; l++)
    {
      const pssalib::datamodel::DataModel_PDM::PropensityIndex & propIdx = ptrPSSACRData->arU3(adjusted_index, l);

      PSSA_TRACE(ptrSimInfo,  << "updating U3 index " << l << " propensity index = (" << propIdx.i << "," << propIdx.j << ")" << std::endl);

      if(adjusted_index == propIdx.i)
        bUpdateSelf = false;

      // Pi
      pssalib::datamodel::CompositionRejectionSamplerData & crsdPi = PSSACRSubVol.crsdPi(propIdx.i);
      const REAL dPi = PSSACRSubVol.arPi(propIdx.i,propIdx.j);
      PSSA_TRACE(ptrSimInfo,  << "dPi = " << dPi << ", min value = " << crsdPi.minValue << std::endl);
      crsdPi.updateValue((UINTEGER)floor(fabs(LOG2(dPi / crsdPi.minValue))) + 1, propIdx.j, dPi);

      // Sigma
      const REAL dSigma = PSSACRSubVol.sigma(propIdx.i);
      PSSA_TRACE(ptrSimInfo,  << "dSigma = " << dSigma << ", inv min sigma = " << invMinSigma << std::endl);
      crsdSigma.updateValue(
        ((dSigma > 0.0) ? ((UINTEGER)floor(fabs(LOG2(dSigma * invMinSigma))) + 1) : 0), propIdx.i, dSigma);
    }

    if(bUpdateSelf)
    {
      const REAL dSigma = PSSACRSubVol.sigma(adjusted_index);
      PSSA_TRACE(ptrSimInfo,  << "self-update for index " << adjusted_index << ": dSigma = " << dSigma << ", inv min sigma = " << invMinSigma << std::endl);
      crsdSigma.updateValue(
          ((dSigma > 0.0) ? ((UINTEGER)floor(fabs(LOG2(dSigma * invMinSigma))) + 1) : 0), adjusted_index, dSigma);
    }
  }

  bool UpdateModule_PSSACR::updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    if(!UpdateModule_PDM::updateSpeciesStructuresReaction(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_PSSACR* ptrPSSACRData =
      static_cast<pssalib::datamodel::DataModel_PSSACR * >
        (ptrSimInfo->getDataModel());

    for(UINTEGER sri = m_sriBegin; sri < m_sriEnd; ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr =
        m_ptrReactionWrapper->getSpeciesReferenceAt(sri);
      if(sr->isReservoir()) continue;
      updateSpeciesStructures(ptrSimInfo, ptrPSSACRData, static_cast<pssalib::datamodel::detail::Subvolume_PSSACR &>(*m_ptrSubvolumeSrc), sr->getIndex());
    }

    return true;
  }

  bool UpdateModule_PSSACR::updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    if(!UpdateModule_PDM::updateSpeciesStructuresDiffusion(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_PSSACR* ptrPSSACRData = 
      static_cast<pssalib::datamodel::DataModel_PSSACR * >
        (ptrSimInfo->getDataModel());

    // update population
    const UINTEGER index = m_ptrReactionWrapper->getSpecies()->getIndex();
    updateSpeciesStructures(ptrSimInfo, ptrPSSACRData,
    static_cast<pssalib::datamodel::detail::Subvolume_PSSACR &>(*m_ptrSubvolumeSrc), index);
    updateSpeciesStructures(ptrSimInfo, ptrPSSACRData,
    static_cast<pssalib::datamodel::detail::Subvolume_PSSACR &>(*m_ptrSubvolumeDst), index);
    return true;
  }

}  } // close namespaces pssalib and update
