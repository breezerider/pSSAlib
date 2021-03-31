/**
 * @file GroupingModule_PDM.cpp
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
 * Grouping module implementation for the Partial-Propensity 
 * Direct Method (Ramaswamy, 2009)
 */

#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/grouping/GroupingModule_PDM.h"

#include "../../include/util/Combinations.h"

namespace pssalib
{
namespace grouping
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  GroupingModule_PDM::GroupingModule_PDM()
  {
    // Do nothing
  }

  //! Copy constructor
  GroupingModule_PDM::GroupingModule_PDM(GroupingModule &g)
    : GroupingModule(g)
  {
    // Do nothing
  }

  //! Destructor
  GroupingModule_PDM::~GroupingModule_PDM()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Initialize data structures (called before each trial)
  bool GroupingModule_PDM::initialize(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    // Call the baseclass method
    if(!GroupingModule::initialize(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_PDM* ptrPDMData = 
      static_cast<pssalib::datamodel::DataModel_PDM*>
        (ptrSimInfo->getDataModel());

    // Check for unsupported reactions
    {
      std::vector<pssalib::datamodel::detail::ReactionWrapper *> arTooManySpecies, arBothReactantsGT1;

      for(UINTEGER rwi = 0; rwi < ptrPDMData->getReactionWrappersCount(); ++rwi)
      {
        pssalib::datamodel::detail::ReactionWrapper & rw = ptrPDMData->getReactionWrapper(rwi);
        if(rw.isDiffusive()) continue; // skip diffusion
        if(rw.getReactantsCount() > 2)
          arTooManySpecies.push_back(&rw);
        else if(rw.getReactantsCount() == 2)
        {
          if((rw.getReactantsListAt(0)->getStoichiometryAbs() > 1)&&
             (rw.getReactantsListAt(1)->getStoichiometryAbs() > 1))
            arBothReactantsGT1.push_back(&rw);
        }
      }

      if((arTooManySpecies.size() > 0) || (arBothReactantsGT1.size() > 0))
      {
        STRINGSTREAM ssTemp;

        if(arTooManySpecies.size() > 0)
        {
          ssTemp << " more than two reactants per reaction "
            "are not supported. Offending reaction are: ";
          std::vector<pssalib::datamodel::detail::ReactionWrapper *>::iterator it = arTooManySpecies.begin();
          for(; it != arTooManySpecies.end(); ++it)
          {
            (*it)->getSymbolicRepresentation(ssTemp);
            ssTemp << "\n";
          }
          PSSA_ERROR(ptrSimInfo, << ssTemp.rdbuf() << std::endl);
        }

        if(arBothReactantsGT1.size() > 0)
        {
          ssTemp << " stoicichiometry of at least one reactant in "
            "every bimolecular reaction must be one. Offending reaction are: ";
          std::vector<pssalib::datamodel::detail::ReactionWrapper *>::iterator it = arBothReactantsGT1.begin();
          for(; it != arBothReactantsGT1.end(); ++it)
          {
            (*it)->getSymbolicRepresentation(ssTemp);
            ssTemp << "\n";
          }
          PSSA_ERROR(ptrSimInfo, << ssTemp.rdbuf() << std::endl);
        }

        bDataLoaded = false;
        return false;
      }
    }

    // Preallocate memory
    UINTEGER l = ptrPDMData->getReactionsCount() / ptrPDMData->getSpeciesCount();
    if(0 == l) l = 1;
    ptrPDMData->arU3.reserve(ptrPDMData->getSpeciesCount() + 1, l);
    ptrPDMData->aruL.reserve(ptrPDMData->getSpeciesCount() + 1, l);

    pssalib::datamodel::detail::JaggedMatrix<UINTEGER> aruLL;
    aruLL.reserve(ptrPDMData->getSpeciesCount() + 1, l);

    ///////////////////////////////////////////////////////
    // Initialize PDM data structures
    pssalib::datamodel::DataModel_PDM::PropensityIndex idxPi;
    boost::scoped_array<REAL> partialPropensity(
      new REAL[ptrPDMData->getSubvolumesCount()]);

    // fill in the mapping data structures
    for(UINTEGER rwi = 0; rwi < ptrPDMData->getReactionWrappersCount(); ++rwi)
    {
      pssalib::datamodel::detail::ReactionWrapper & rw = ptrPDMData->getReactionWrapper(rwi);

      // store the rate
      idxPi.rate = ptrPDMData->getReactionWrapper(rwi).getRate();

      PSSA_TRACE(ptrSimInfo, << "= reaction : " << rw.toString() << std::endl);

      // specific probability rate
      std::fill_n(partialPropensity.get(), ptrPDMData->getSubvolumesCount(), rw.getRate());

      bool selfDep = false;
      if(rw.isDiffusive())
      {
        idxPi.i = rw.getSpecies()->getIndex() + 1;

        // compute partial propensities
        for(UINTEGER svi = 0; svi < ptrPDMData->getSubvolumesCount(); ++svi)
          partialPropensity[svi] *= 2.0 * (REAL)ptrPDMData->getDimsCount();

//         PSSA_TRACE(ptrSimInfo, << "== diffusion assigned to species #" << idxPi.i
//           << " : D=" << rw.getSpecies()->getDiffusionConstant() << std::endl);
      }
      else
      {
        if(rw.getReactantsCount() > 1) // affects exactly two species (see above)
        {
          const pssalib::datamodel::detail::SpeciesReference * sr1 = rw.getReactantsListAt(0);
          const pssalib::datamodel::detail::SpeciesReference * sr2 = rw.getReactantsListAt(1);

          if(1 != sr2->getStoichiometryAbs())
          {
            rw.swapSpeciesReferencesAt(0,1);
            sr1 = rw.getReactantsListAt(0);
            sr2 = rw.getReactantsListAt(1);
            PSSA_TRACE(ptrSimInfo, << "== swapping species in reaction #" << rwi << std::endl);
          }

          // compute partial propensities
          for(UINTEGER svi = 0; svi < ptrPDMData->getSubvolumesCount(); ++svi)
            partialPropensity[svi] *= pssalib::util::getPartialCombinationsHeteroreactions(ptrPDMData->getSubvolume(svi).population(sr1->getIndex()), sr1->getStoichiometryAbs());

          // row index in PI for this reaction
          idxPi.i = sr2->getIndex() + 1;
          // column index in PI for this reaction
          idxPi.j = ptrPDMData->aruL.get_cols(idxPi.i);
          // dependent species stoichiometry
          idxPi.stoichiometry = sr1->getStoichiometryAbs();

          // store dependency
          ptrPDMData->arU3.push_back(sr1->getIndex() + 1, idxPi);
        }
        else // involves a single species
        {
          const pssalib::datamodel::detail::SpeciesReference * sr1 = rw.getReactantsListAt(0);

          if(sr1->isReservoir())
            idxPi.i = 0;
          else
          {
            // row index in PI for this reaction
            idxPi.i = sr1->getIndex() + 1;

            // compute partial propensities
            for(UINTEGER svi = 0; svi < ptrPDMData->getSubvolumesCount(); ++svi)
              partialPropensity[svi] *= pssalib::util::getPartialCombinationsHomoreactions(ptrPDMData->getSubvolume(svi).population(sr1->getIndex()), sr1->getStoichiometryAbs());

            // Account for self dependency
            if(sr1->getStoichiometryAbs() > 1)
            {
              selfDep = true;
              // column index in PI for this reaction
              idxPi.j = ptrPDMData->aruL.get_cols(idxPi.i);
              // dependent species stoichiometry
              idxPi.stoichiometry = sr1->getStoichiometryAbs();

              ptrPDMData->arU3.push_back(idxPi.i, idxPi);
            }
          }
        }
      }

      PSSA_TRACE(ptrSimInfo, << "== " << ((rw.getReactantsCount() > 0) ? (((selfDep)||(rw.getReactantsCount() > 1)) ? "bimolecular" : "unimolecular") : "diffusion") << ((idxPi.i > 0) ? " assigned to species #" : " assigned to reservoir species ") << ((idxPi.i > 0) ? idxPi.i-1 : 0) << " with partial propensity = " << partialPropensity[0] << std::endl);

      // position in PI --> reaction number
      pssalib::datamodel::detail::ReactionWrapper * ptrRW = &ptrPDMData->getReactionWrapper(rwi);
      ptrPDMData->aruL.push_back(idxPi.i, ptrRW);
      UINTEGER rwi1 = rwi + 1;
      aruLL.push_back(idxPi.i, rwi1);

      // assign partial propensities
      for(UINTEGER svi = 0; svi < ptrPDMData->getSubvolumesCount(); ++svi)
      {
        ptrPDMData->getSubvolume(svi).arPi.push_back(idxPi.i, partialPropensity[svi]);
        ptrPDMData->getSubvolume(svi).lambda(idxPi.i) += partialPropensity[svi];
      }
    }

    PSSA_TRACE(ptrSimInfo, << "Mapping variables ready.\naruL : " << aruLL << "\narU3 : \n" << ptrPDMData->arU3 << "\nSample arPi : \n" << ptrPDMData->getSubvolume(0).arPi << std::endl);

    ptrPDMData->dTotalPropensity = 0.0;
    for(UINTEGER svi = 0; svi < ptrPDMData->getSubvolumesCount(); ++svi)
    {
      ptrPDMData->getSubvolume(svi).dTotalPropensity = 0.0;
      for(UINTEGER si = 0; si < ptrPDMData->getSpeciesCount() + 1; ++si)
      {
        if(0 == si) // reservoir species
        {
          ptrPDMData->getSubvolume(svi).sigma(si) = ptrPDMData->getSubvolume(svi).lambda(si);

          PSSA_TRACE(ptrSimInfo, << "== Reservoir species : Lambda  = "
            << ptrPDMData->getSubvolume(svi).lambda(si) << "; Sigma = "
            << ptrPDMData->getSubvolume(svi).sigma(si) << std::endl);
        }
        else
        {
          ptrPDMData->getSubvolume(svi).sigma(si) = ptrPDMData->getSubvolume(svi).population(si-1) * ptrPDMData->getSubvolume(svi).lambda(si);

          PSSA_TRACE(ptrSimInfo, << "== Species #" << si-1 << " '" 
            << ptrPDMData->getSpecies(si-1)->toString() << "' : Lambda  = "
            << ptrPDMData->getSubvolume(svi).lambda(si) << "; Sigma = "
            << ptrPDMData->getSubvolume(svi).sigma(si) << std::endl);
        }

        ptrPDMData->getSubvolume(svi).dTotalPropensity += ptrPDMData->getSubvolume(svi).sigma(si);
      }
      PSSA_TRACE(ptrSimInfo, << "= total propensity = " 
        << ptrPDMData->getSubvolume(svi).dTotalPropensity << std::endl);
      ptrPDMData->dTotalPropensity += ptrPDMData->getSubvolume(svi).dTotalPropensity;
    }

    PSSA_TRACE(ptrSimInfo, << "global total propensity = " << ptrPDMData->dTotalPropensity << std::endl);

    return true;
  }

}  } // close namespaces pssalib and grouping
