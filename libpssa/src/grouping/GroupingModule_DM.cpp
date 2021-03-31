/**
 * @file GroupingModule_DM.cpp
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
 * Grouping module implementation for the Gillespie's Direct Method
 */

#include "../../include/datamodel/DataModel_DM.h"
#include "../../include/grouping/GroupingModule_DM.h"

#include "../../include/util/Combinations.h"

namespace pssalib
{
namespace grouping
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  GroupingModule_DM::GroupingModule_DM()
  {
    // Do nothing
  }

  //! Copy constructor
  GroupingModule_DM::GroupingModule_DM(GroupingModule &g)
    : GroupingModule(g)
  {
    // Do nothing
  }

  //! Destructor
  GroupingModule_DM::~GroupingModule_DM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Initialize data structures (called before each trial)
  bool GroupingModule_DM::initialize(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    // Call the baseclass method
    if(!GroupingModule::initialize(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_DM* ptrDMData = 
      static_cast<pssalib::datamodel::DataModel_DM * >
        (ptrSimInfo->getDataModel());

    for(UINTEGER svi = 0; svi < ptrDMData->getSubvolumesCount(); ++svi)
    {
      pssalib::datamodel::detail::Subvolume_DM & DMSubVol = ptrDMData->getSubvolume(svi);

      DMSubVol.dTotalPropensity = 0.0;

      // Fill the propensities array
      // Normal reactions
      for(UINTEGER rwi = 0; rwi < ptrDMData->getReactionWrappersCount(); rwi++)
      {
        pssalib::datamodel::detail::ReactionWrapper & rw = ptrDMData->getReactionWrapper(rwi);

        DMSubVol.propensity(rwi) = 0.0;

        // compute reaction propensity
        REAL temp = rw.getRate();
        if(rw.isDiffusive())
        {
          const pssalib::datamodel::detail::Species * sp = rw.getSpecies();
          temp *= (REAL)DMSubVol.population(sp->getIndex()) * 2.0 * (REAL)ptrDMData->getDimsCount();
        }
        else
        {
          for(UINTEGER ri = 0; ri < rw.getReactantsCount(); ++ri)
          {
            const pssalib::datamodel::detail::SpeciesReference * sr = rw.getReactantsListAt(ri);
            if(!sr->isReservoir())
              temp *= pssalib::util::getPartialCombinationsHeteroreactions(DMSubVol.population(sr->getIndex()), sr->getStoichiometryAbs());
          }
        }

        // store propensity on the subvolume scale
        DMSubVol.propensity(rwi) = temp;
        PSSA_TRACE(ptrSimInfo, << "propensity_" << rwi << " = " << temp  << "; from array = " << DMSubVol.propensity(rwi) << std::endl);
        DMSubVol.dTotalPropensity += temp;
      }

      PSSA_TRACE(ptrSimInfo, << "totalPropensity=" << DMSubVol.dTotalPropensity << std::endl);

      // update global structures
      ptrDMData->dTotalPropensity += DMSubVol.dTotalPropensity;
    }

    return true;
  }
}
}
