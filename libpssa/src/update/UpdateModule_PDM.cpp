/**
 * @file UpdateModule_PDM.cpp
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
 * Update module implementation for the Partial-Propensity 
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/detail/SpeciesReference.h"
#include "../../include/update/UpdateModule_PDM.h"

#include "../../include/util/Combinations.h"

namespace pssalib
{
namespace update
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  UpdateModule_PDM::UpdateModule_PDM()
  {
    // Do nothing
  }

  //! Destructor
  UpdateModule_PDM::~UpdateModule_PDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Update reaction propensities
  bool UpdateModule_PDM::updateSpeciesStructures(pssalib::datamodel::SimulationInfo * ptrSimInfo,
                                                 pssalib::datamodel::DataModel_PDM * ptrPDMData,
                                                 pssalib::datamodel::detail::Subvolume_PDM & PDMSubVol,
                                                 UINTEGER index, INTEGER stoichiometry)
  {
    UINTEGER population = PDMSubVol.population(index);
    REAL dTotalPropensityChange = 0.0;
    bool updateSelf = true;

    // For each reaction that contains current speciesIndex
    for(UINTEGER l = 0, U3_rowlen = ptrPDMData->arU3.get_cols(index+1); l < U3_rowlen; l++)
    {
      const pssalib::datamodel::DataModel_PDM::PropensityIndex & propIdx =
        ptrPDMData->arU3(index+1, l);

      REAL newProp = propIdx.rate, oldProp = PDMSubVol.arPi(propIdx.i,propIdx.j), temp;

      // take polymerization into account
      if(propIdx.i != (index+1))
      {
        newProp *= ((REAL)pssalib::util::getPartialCombinationsHeteroreactions(population, propIdx.stoichiometry));

        temp = newProp - oldProp;

        PDMSubVol.arPi(propIdx.i,propIdx.j) = newProp;
        PDMSubVol.lambda(propIdx.i) += temp;

        temp = PDMSubVol.population(propIdx.i-1) * PDMSubVol.lambda(propIdx.i);

        dTotalPropensityChange += temp - PDMSubVol.sigma(propIdx.i);
        PDMSubVol.sigma(propIdx.i) = temp;
      }
      else
      {
        updateSelf = false;
        newProp *= ((REAL)pssalib::util::getPartialCombinationsHomoreactions(population, propIdx.stoichiometry));

        temp = newProp - oldProp;

        PDMSubVol.arPi(propIdx.i,propIdx.j) = newProp;
        PDMSubVol.lambda(propIdx.i) += temp;

        temp = population * PDMSubVol.lambda(propIdx.i);

        dTotalPropensityChange += temp - PDMSubVol.sigma(propIdx.i);
        PDMSubVol.sigma(propIdx.i) = temp;
      }

      PSSA_TRACE(ptrSimInfo, << "updating reaction '" << ptrPDMData->aruL(propIdx.i, propIdx.j)->toString() << "' affected by species #"
        << index << " (prop index=" << propIdx.i << "; stoichiometry=" << propIdx.stoichiometry << "; population=" << population << "; change=" << stoichiometry
        << ") with old pp = " << oldProp << " and new pp = " << newProp << (updateSelf ? "; self-update required" : "") << std::endl);

      PSSA_TRACE(ptrSimInfo, << "lambda[" << propIdx.i << "]" << PDMSubVol.lambda(propIdx.i) << "\tsigma[" << propIdx.i << "]" << PDMSubVol.sigma(propIdx.i) << "\tarPi :\n" << PDMSubVol.arPi << std::endl);
    }

    // Update the group propensity of the affected species
    // (if it was not updated already)
    if( updateSelf )
    {
      REAL temp = population * PDMSubVol.lambda(index+1);
      dTotalPropensityChange += temp - PDMSubVol.sigma(index+1);
      PDMSubVol.sigma(index+1) = temp;
    }

    PDMSubVol.dTotalPropensity += dTotalPropensityChange;

    return true;
  }

  bool UpdateModule_PDM::updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel_PDM * ptrPDMData =
      static_cast<pssalib::datamodel::DataModel_PDM * >
        (ptrSimInfo->getDataModel());

    for(UINTEGER sri = m_sriBegin; sri < m_sriEnd; ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * sr =
        m_ptrReactionWrapper->getSpeciesReferenceAt(sri);
      if(sr->isConstant()) continue;
      if(!updateSpeciesStructures(ptrSimInfo, ptrPDMData,
        static_cast<pssalib::datamodel::detail::Subvolume_PDM &>(*m_ptrSubvolumeSrc),
        sr->getIndex(), ((sri < m_sriReactants) ? -1 : 1) * sr->getStoichiometry()))
        return false;
    }

    return true;
  }

  bool UpdateModule_PDM::updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel_PDM * ptrPDMData = 
      static_cast<pssalib::datamodel::DataModel_PDM * >
        (ptrSimInfo->getDataModel());

    // update population
    const UINTEGER index = m_ptrReactionWrapper->getSpecies()->getIndex();
    return (updateSpeciesStructures(ptrSimInfo, ptrPDMData,
      static_cast<pssalib::datamodel::detail::Subvolume_PDM &>(*m_ptrSubvolumeSrc), index, -1) &&
      updateSpeciesStructures(ptrSimInfo, ptrPDMData,
      static_cast<pssalib::datamodel::detail::Subvolume_PDM &>(*m_ptrSubvolumeDst), index, 1));
  }

}  } // close namespaces pssalib and update
