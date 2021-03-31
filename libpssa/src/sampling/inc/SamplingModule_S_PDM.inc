/**
 * @file SamplingModule_S_PDM.inc
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
 * Implementation of the sampling for the Partial-Propensity 
 * Direct Method and the Sorting Partial Propensity Direct
 * Method (Ramaswamy, 2009)
 */


#ifdef PSSALIB_INTERNAL_S_PDM_MAP_I
#undef PSSALIB_INTERNAL_S_PDM_MAP_I
#endif
#ifdef PSSALIB_INTERNAL_S_PDM_MAP_IJ
#undef PSSALIB_INTERNAL_S_PDM_MAP_IJ
#endif

#if defined(PSSALIB_INTERNAL_PDM_MODULE)
#  define PSSALIB_INTERNAL_S_PDM_MAP_I(vol,i) i
#  define PSSALIB_INTERNAL_S_PDM_MAP_IJ(vol,i,j) j
#elif defined(PSSALIB_INTERNAL_SPDM_MODULE)
#  define PSSALIB_INTERNAL_S_PDM_MAP_I(vol,i) vol.mapRowIndex(i)
#  define PSSALIB_INTERNAL_S_PDM_MAP_IJ(vol,i,j) vol.mapColIndex(i, j)
#endif

////////////////////////////////
// PDM & SPDM sampling
{
  // Cast the data model to a suitable type
#if defined(PSSALIB_INTERNAL_PDM_MODULE)
  pssalib::datamodel::DataModel_PDM * 
#elif defined(PSSALIB_INTERNAL_SPDM_MODULE)
  pssalib::datamodel::DataModel_SPDM * 
#endif
    ptrData =
#if defined(PSSALIB_INTERNAL_PDM_MODULE)
    static_cast<pssalib::datamodel::DataModel_PDM *>(ptrSimInfo->getDataModel())
#elif defined(PSSALIB_INTERNAL_SPDM_MODULE)
    static_cast<pssalib::datamodel::DataModel_SPDM *>(ptrSimInfo->getDataModel())
#endif
    ;

#if defined(PSSALIB_INTERNAL_PDM_MODULE)
  pssalib::datamodel::detail::Subvolume_PDM &
#elif defined(PSSALIB_INTERNAL_SPDM_MODULE)
  pssalib::datamodel::detail::Subvolume_SPDM &
#endif
    SubVol = ptrData->getSubvolume(ptrData->nu);

  // Temporary variables
  REAL temp1, temp2;
  UINTEGER i, j, N = ptrData->getSpeciesCount() + 1, tempI, tempJ;

  // Sample reaction
  temp1 = gsl_rng_uniform_pos (m_ptrRNG) * SubVol.dTotalPropensity;
  temp2 = 0.0;
  for(i = 0; i < N; i++)
  {
    tempI = PSSALIB_INTERNAL_S_PDM_MAP_I(SubVol,i);
    temp2 += SubVol.sigma(tempI);

    if(temp1 < temp2)
      break;
  }

  if( i >= N )
  {
    if (temp1/SubVol.dTotalPropensity > 1.01) // exclude numeric noise
    {
      PSSA_ERROR(ptrSimInfo, << "row target is above cumulative sum of propensities: "
        << temp1 << " > " << temp2 << std::endl);
      return false;
    }
    else
    {
      do
      {
        i--; tempI = PSSALIB_INTERNAL_S_PDM_MAP_I(SubVol,i);
      } while((i > 0)&&(0.0 >= SubVol.sigma(tempI)));
    }
    PSSA_WARNING(ptrSimInfo, << "row target is too close to cumulative sum of propensities: "
      << temp1 << " ~ " << temp2 << "; choosing entry #" << i << " ('" 
      << ptrData->getSpecies(i)->toString() << "') as fail-safe "
      "(sigma " << SubVol.sigma(tempI) << ")" << std::endl);
  }
// #ifdef DEBUG
//   PSSA_TRACE(ptrSimInfo, << "row target=" << temp1 << "; row target max=" << temp2 << "; tempI=" << tempI << "; population=" << ((0 != tempI) ? SubVol.population(tempI-1) : 1) << std::endl);
// #endif
  temp1 = (temp1 - temp2 + SubVol.sigma(tempI));
  if(0 != tempI)
    temp1 /= (REAL)SubVol.population(tempI-1);
  temp2 = 0.0;
  N = SubVol.arPi.get_cols(tempI);
  for(j = 0; j < N; j++)
  {
    tempJ = PSSALIB_INTERNAL_S_PDM_MAP_IJ(SubVol,i,j);
    temp2 += SubVol.arPi(tempI, tempJ);

    if(temp1 < temp2)
      break;
  }

  if( j >= N )
  {
    if (temp1/SubVol.dTotalPropensity > 1.01) // exclude numeric noise
    {
      PSSA_ERROR(ptrSimInfo, << "column target is above cumulative sum of propensities: "
        << temp1 << " > " << temp2 << std::endl);
      return false;
    }
    else
    {
      do
      {
        j--; tempJ = PSSALIB_INTERNAL_S_PDM_MAP_IJ(SubVol,i,j);
      } while((j > 0)&&(0.0 >= SubVol.arPi(tempI, tempJ)));
      PSSA_WARNING(ptrSimInfo, << "col target is too close to cumulative sum of propensities: "
        << temp1 << " ~ " << temp2 << "; choosing entry #" << j << " as a fail-safe." << std::endl);
    }
  }

#if defined(PSSALIB_INTERNAL_SPDM_MODULE)
  ptrData->rowIndex = i;
  ptrData->colIndex = j;
#endif

  // Set the next reaction index
  ptrData->mu = ptrData->aruL(tempI, tempJ)->getSerialNumber();

  return true;
}

#ifdef PSSALIB_INTERNAL_S_PDM_MAP_I
#undef PSSALIB_INTERNAL_S_PDM_MAP_I
#endif
#ifdef PSSALIB_INTERNAL_S_PDM_MAP_IJ
#undef PSSALIB_INTERNAL_S_PDM_MAP_IJ
#endif
