/**
 * @file UpdateModule_DM.cpp
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
 * Update module implementation for the Gillespie's Direct Method
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/update/UpdateModule_DM.h"

#include "../../include/util/Combinations.h"

namespace pssalib
{
namespace update
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  UpdateModule_DM::UpdateModule_DM()
  {
    // Do nothing
  }

  //! Destructor
  UpdateModule_DM::~UpdateModule_DM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Update reaction propensities
  bool UpdateModule_DM::updateSpeciesStructures(pssalib::datamodel::SimulationInfo * ptrSimInfo,
                                                pssalib::datamodel::DataModel_DM * ptrDMData,
                                                pssalib::datamodel::detail::Subvolume_DM & DMSubVol)
  {
    DMSubVol.dTotalPropensity = 0.0;

    // Recalculate the propensity array
    for(UINTEGER rwi = 0; rwi < ptrDMData->getReactionWrappersCount(); ++rwi)
    {
      pssalib::datamodel::detail::ReactionWrapper & rw =
        ptrDMData->getReactionWrapper(rwi);

      DMSubVol.propensity(rwi) = 0.0;

      PSSA_TRACE(ptrSimInfo,  << "updating reaction index " << rwi << std::endl);

      // compute reaction propensity
      REAL temp = rw.getRate();
      if(rw.isDiffusive())
      {
        const pssalib::datamodel::detail::Species * sp = rw.getSpecies();
        temp *= DMSubVol.population(sp->getIndex()) * 2.0 * (REAL)ptrDMData->getDimsCount();
      }
      else
      {
        for(UINTEGER ri = 0; ri < rw.getReactantsCount(); ++ri)
        {
          const pssalib::datamodel::detail::SpeciesReference * sr = rw.getReactantsListAt(ri);
          PSSA_TRACE(ptrSimInfo, << "updating species index " << sr->getIndex() << " with stoichiometry " << sr->getStoichiometryAbs() << std::endl);
          if(!sr->isReservoir())
            temp *= pssalib::util::getPartialCombinationsHeteroreactions(DMSubVol.population(sr->getIndex()), sr->getStoichiometryAbs());
        }
      }

      // store propensity on the subvolume scale
      DMSubVol.propensity(rwi) = temp;
      DMSubVol.dTotalPropensity += temp;
    }

    return true;
  }

  bool UpdateModule_DM::updateSpeciesStructuresReaction(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel_DM* ptrDMData =
      static_cast<pssalib::datamodel::DataModel_DM * >
        (ptrSimInfo->getDataModel());

    return updateSpeciesStructures(ptrSimInfo, ptrDMData, ptrDMData->getSubvolume(ptrDMData->nu));
  }

  bool UpdateModule_DM::updateSpeciesStructuresDiffusion(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel_DM* ptrDMData = 
      static_cast<pssalib::datamodel::DataModel_DM * >
        (ptrSimInfo->getDataModel());
    return updateSpeciesStructures(ptrSimInfo, ptrDMData, ptrDMData->getSubvolume(ptrDMData->nu)) &&
      updateSpeciesStructures(ptrSimInfo, ptrDMData, ptrDMData->getSubvolume(ptrDMData->nu_D));
  }

}  } // close namespaces pssalib and update
