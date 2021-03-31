/**
 * @file GroupingModule_PSSACR.cpp
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
 * Grouping module implementation for the Partial Propensity SSA 
 * with Composition-Rejection sampling (Ramaswamy, 2010)
 */

#include "../../include/datamodel/DataModel_PSSACR.h"
#include "../../include/datamodel/PSSACR_Bins.h"
#include "../../include/grouping/GroupingModule_PSSACR.h"

#include "../../include/util/Combinations.h"
#include "../../include/util/Maths.h"

namespace pssalib
{
namespace grouping
{
  ////////////////////////////////
  // Constructors

  //! Default constructor
  GroupingModule_PSSACR::GroupingModule_PSSACR()
  {
    // Do nothing
  }

  //! Copy constructor
  GroupingModule_PSSACR::GroupingModule_PSSACR(GroupingModule & g)
    : GroupingModule_PDM(g)
  {
    // Do nothing
  }

  //! Destructor
  GroupingModule_PSSACR::~GroupingModule_PSSACR()
  {
    // Do nothing
  }

  ////////////////////////////////
  // Methods

  //! Initialize data structures (called before each trial)
  bool GroupingModule_PSSACR::initialize(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    // Call the base class method
    if(!GroupingModule_PDM::initialize(ptrSimInfo))
      return false;

    // Cast the data model to a suitable type
    pssalib::datamodel::DataModel_PSSACR* ptrPSRDCRData = 
      static_cast<pssalib::datamodel::DataModel_PSSACR *>
        (ptrSimInfo->getDataModel());

    // Temporary variables
    REAL minSigma = std::numeric_limits<REAL>::max();
    boost::scoped_array<REAL> minPi((new REAL[ptrPSRDCRData->getSpeciesCount()]));

    ////////////////////////////////////////////////
    // Calculate the minimum values
    bool bSetSigma = false;
    for(UINTEGER si = 0, unPi; si < ptrPSRDCRData->getSpeciesCount(); ++si)
    {
      unPi = ptrPSRDCRData->aruL.get_cols(si);

      if(0 == unPi) // not a reactant
        minPi[si] = 0.0;
      else
        minPi[si] = std::numeric_limits<REAL>::max();

      for(UINTEGER sj = 0; sj < unPi; ++sj)
      {
        pssalib::datamodel::detail::ReactionWrapper * rw = ptrPSRDCRData->aruL(si,sj);
        
        REAL temp = rw->getRate();
        const pssalib::datamodel::detail::SpeciesReference * sr =
          rw->getReactantsListAt(0);
        if(sr->getIndex() == si)
        {
          temp *= pssalib::util::getPartialCombinationsHomoreactions(sr->getStoichiometryAbs(), sr->getStoichiometryAbs());
        }
        else
        {
          temp *= pssalib::util::getPartialCombinationsHeteroreactions(sr->getStoichiometryAbs(), sr->getStoichiometryAbs());
        }

        if(minPi[si] > temp)
          minPi[si] = temp;

        if(temp > 0.0)
        {
          bSetSigma = true;

          if((sr->getIndex() == si)&&(sr->getStoichiometryAbs() > 0))
            temp *= (REAL)sr->getStoichiometryAbs();

          if(minSigma > temp)
            minSigma = temp;
        }
      }
    }
    
    if(!bSetSigma) minSigma = 0.0;

    ////////////////////////////////////////////////
    // Compute distribution
    for(UINTEGER svi = 0; svi < ptrPSRDCRData->getSubvolumesCount(); ++svi)
    {
      pssalib::datamodel::detail::Subvolume_PSSACR & PSSACRSubVol = ptrPSRDCRData->getSubvolume(svi);
      
      PSSACRSubVol.crsdSigma.minValue = minSigma;
      PSSACRSubVol.crsdSigma.bins.clear();
      PSSACRSubVol.crsdSigma.bins.resize(ptrPSRDCRData->getSpeciesCount());
      for(UINTEGER si = 0, unPi, k; si < ptrPSRDCRData->getSpeciesCount(); ++si)
      {
        PSSACRSubVol.crsdPi(si).minValue = minPi[si];
        
        // Sigma's
        if(0 != PSSACRSubVol.sigma(si))
        {
          if (ptrPSRDCRData->crsdVolume.minValue == 0.0)
            k = (UINTEGER)floor(fabs(LOG2(PSSACRSubVol.sigma(si)))) + 1;
          else
            k = (UINTEGER)floor(fabs(LOG2(PSSACRSubVol.sigma(si) / PSSACRSubVol.crsdSigma.minValue))) + 1;
          PSSACRSubVol.crsdSigma.bins.updateValue(k, si, PSSACRSubVol.sigma(si));
        }

        // Pi's
        unPi = PSSACRSubVol.arPi.get_cols(si);
        PSSACRSubVol.crsdPi(si).bins.resize(unPi);
        for(UINTEGER sj = 0; sj < unPi; ++sj)
        {
          if(0 == PSSACRSubVol.arPi(si,sj))
            continue;

          k = (UINTEGER)floor(fabs(LOG2(PSSACRSubVol.arPi(si,sj) / PSSACRSubVol.crsdPi(si).minValue))) + 1;

          //k = (UINTEGER)std::floor(log2(ptrPSRDCRData->arSubvolumes_PDM[si].arPi(i,j) / ptrPSRDCRData->arSubvolumes_PSSACR[si].crsdPi[i].minValue)) + 1;
          PSSACRSubVol.crsdPi(si).bins.updateValue(k, sj, PSSACRSubVol.arPi(si,sj));
        }
      }
    }

    return true;
  }

}  } // close namespaces pssalib and grouping
