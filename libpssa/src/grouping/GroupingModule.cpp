/**
 * @file GroupingModule.cpp
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
 * Implementation of the generic Grouping module:
 *   - validate and parse SBML model
 *   - fill in method-specific simulation data structures with information from
 *     the SBML model
 *   - prepare the datastructures before a simulation starts
 */

#include "../../include/datamodel/SimulationInfo.h"
#include "../../include/grouping/GroupingModule.h"

#include "../../include/util/Indexing.h"
#include "../../include/util/Maths.h"
#include "../../include/util/IO.hpp"

#include <boost/iostreams/filtering_stream.hpp>

namespace pssalib
{
namespace grouping
{
  ////////////////////////////////
  // Constructors

  //! Constructor
  GroupingModule::GroupingModule()
  {
    // Do nothing
  }

  //! Copy constructor
  GroupingModule::GroupingModule(GroupingModule & other)
    : bDataLoaded(other.bDataLoaded)
  {
    // Do nothing
  }

  //! Destructor
  GroupingModule::~GroupingModule()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  bool GroupingModule::preinitialize(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel* ptrData = ptrSimInfo->getDataModel();

    if(NULL == ptrData)
    {
      PSSA_ERROR(ptrSimInfo, << "data structures not initialized!\n");
      return false;
    }

    // Clean-up previous model
    ptrData->clear();

    bool bResult = true;
    try
    {
      // make a shallow copy of the user defined model
      ((pssalib::datamodel::detail::Model*)ptrData)->copy(ptrSimInfo->getModel());

      ptrData->setup(ptrSimInfo->getDimsCount(), ptrSimInfo->getDims(), ptrSimInfo->eBoundaryConditions);
      PSSA_TRACE(ptrSimInfo, << "done setting up volume decomposition.\n");

      ptrData->nu = 0; ptrData->nu_D = 0;

      STRINGSTREAM ssTemp;
      boost::iostreams::filtering_ostream os;
      os.push(pssalib::io::prefix_filter("(INFO) : "));
      os.push(ssTemp);

      ptrData->printReactionNetwork(os);
      os.flush();

      PSSA_INFO(ptrSimInfo, << ssTemp.rdbuf() << std::flush);
    }
    catch(std::runtime_error & e)
    {
      PSSA_ERROR(ptrSimInfo, << "processing model '" << ptrData->getId() << "' failed: " << e.what() << std::endl);
      bResult = false;
    }
    catch(std::bad_alloc & e)
    {
      PSSA_ERROR(ptrSimInfo, << "ran out of memory while processing model '" << ptrData->getId() << "': " << e.what() << std::endl);
      bResult = false;
    }

    // Set the global flag
    bDataLoaded = bResult;
    return bResult;
  }

  //! Initialize data structures (called before each trial)
  bool GroupingModule::initialize(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel * ptrData = ptrSimInfo->getDataModel();

    // clear all data structures
    ptrData->clear();
    ptrData->vQueuedReactions.clear();
    ptrData->dTotalPropensity = 0.0;

    // setup the initial population
    if(ptrData->getSpeciesCount() > 0)
    {
      boost::scoped_array< UINTEGER > arPopulation(new UINTEGER[ptrData->getSubvolumesCount() * ptrData->getSpeciesCount()]);
      boost::scoped_array< UINTEGER * > arPtrPopulation(new UINTEGER *[ptrData->getSubvolumesCount()]);
      for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
        arPtrPopulation[svi] = arPopulation.get() + svi * ptrData->getSpeciesCount();
      memset(arPopulation.get(), 0, sizeof(UINTEGER) * ptrData->getSubvolumesCount() * ptrData->getSpeciesCount());

      switch(ptrSimInfo->eInitialPopulation)
      {
        case pssalib::datamodel::detail::IP_Distribute:
        {
          REAL invSubvolCount = REAL(1.0)/REAL(ptrData->getSubvolumesCount());
          for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
            for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
              arPtrPopulation[svi][si] = std::floor(REAL(ptrData->getSpecies(si)->getInitialAmount()) * invSubvolCount);
        }
        break;
        case pssalib::datamodel::detail::IP_Concentrate:
        {
          boost::scoped_array<UINTEGER> subMid(new UINTEGER[ptrData->getDimsCount()]);
          for(BYTE d = 0; d < ptrData->getDimsCount(); ++d)
            subMid[d] = ptrData->getDims(d) / 2;
          UINTEGER idxMid = 0;
          util::sub2ind(ptrData->getDimsCount(), ptrData->getDims(), subMid.get(), idxMid);

          for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
            arPtrPopulation[idxMid][si] = ptrData->getSpecies(si)->getInitialAmount();
        }
        break;
        case pssalib::datamodel::detail::IP_Multiply:
        {
          for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
            for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
              arPtrPopulation[svi][si] = ptrData->getSpecies(si)->getInitialAmount();
        }
        break;
        case pssalib::datamodel::detail::IP_UserDefined:
        {
          if(NULL != ptrSimInfo->ptrPopulationInitializer)
            (*ptrSimInfo->ptrPopulationInitializer)(ptrData, arPtrPopulation.get(), ptrSimInfo->ptrPopulationInitializerUserData);
          else
          {
            PSSA_ERROR(ptrSimInfo, << "when using a user-defined initial population initializer, the ptrInitialPopulationUserDefined must not be NULL and contain a pointer to the initializer function" << std::endl);
            return false;
          }
        }
        break;
        case pssalib::datamodel::detail::IP_Invalid:
        default:
        {
          if(ptrData->getDimsCount() > 1)
          {
            PSSA_ERROR(ptrSimInfo, << "initial population initializer undefined" << std::endl);
            return false;
          }
          else
          {
            for(UINTEGER si = 0; si < ptrData->getSpeciesCount(); ++si)
              arPtrPopulation[0][si] = ptrData->getSpecies(si)->getInitialAmount();
          }
        }
        break;
      }
      ptrData->setupPopulation(arPtrPopulation.get());
    }

    return true;
  }

  void GroupingModule::postInitialize(pssalib::datamodel::SimulationInfo * ptrSimInfo)
  {
    pssalib::datamodel::DataModel * ptrData = ptrSimInfo->getDataModel();

    ptrData->crsdVolume.bins.resize(ptrData->getSubvolumesCount());
    // Distribute the subvolume propensities into the bins.
    for(UINTEGER svi = 0; svi < ptrData->getSubvolumesCount(); ++svi)
    {
      const pssalib::datamodel::detail::Subvolume & sv = ptrData->getSubvolume(svi);
      UINTEGER k = (UINTEGER)std::floor(fabs(LOG2(sv.dTotalPropensity /
                                                 ptrData->crsdVolume.minValue))) + 1;
      PSSA_TRACE(ptrSimInfo, << "Subvol #" << svi << " : tot_prop="
        << sv.dTotalPropensity << "; k=" << k << std::endl);
      ptrData->crsdVolume.updateValue(k, svi, sv.dTotalPropensity);
    }
  }

}  } // close namespaces pssalib and grouping
