/**
 * @file PSSACR_Bins.cpp
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
 * Implements methods for propensity binning when using PSSA-CR
 */

#include "../../include/stdheaders.h"
#include "../../include/typedefs.h"
#include "../../include/datamodel/PSSACR_Bins.h"

namespace pssalib
{
namespace datamodel
{
  ////////////////////////////////
  // Bin class

  //! Constructor
  PSSACR_Bin::PSSACR_Bin()
    : dBinSum(0.0)
    , arunBinEl(NULL)
    , unNumBinEl(0)
    , unCapBinEl(0)
  {
    // Do nothing
  }

  //! Copy constructor
  PSSACR_Bin::PSSACR_Bin(const PSSACR_Bin &b)
    : dBinSum(0.0)
    , arunBinEl(NULL)
    , unNumBinEl(0)
    , unCapBinEl(0)
  {
    if(NULL != b.arunBinEl)
    {
      resize(b.unCapBinEl);
      unNumBinEl = b.unNumBinEl;
      memcpy(b.arunBinEl, arunBinEl, unNumBinEl*sizeof(UINTEGER));
      dBinSum    = b.dBinSum;
    }
  }

  //! Destructor
  PSSACR_Bin::~PSSACR_Bin()
  {
    // Free memory
    clear();
  }

  //! Resize the bin
  void PSSACR_Bin::resize(UINTEGER unsize)
  {
    clear();
    if(unsize > 0)
    {
      arunBinEl = new UINTEGER[unsize];
      unCapBinEl = unsize;
    }
  }

  //! Clear the bin
  void PSSACR_Bin::clear()
  {
    if(NULL != arunBinEl)
    {
      delete [] arunBinEl;
      arunBinEl  = NULL;
      unCapBinEl = 0;
      unNumBinEl = 0;
      dBinSum    = 0.0;
    }
  }

  /**
   * Add an index to the list
   * 
   * @param el propensity index to be added to this bin
   * @return index of the newly added element
   */
  UINTEGER PSSACR_Bin::push_back(UINTEGER el)
  {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(unNumBinEl >= unCapBinEl)
        throw std::runtime_error("PSSACR_Bin::push_back() - number of elements exceeds vector capacity.");
      else
#endif
        arunBinEl[unNumBinEl] = el;
    return unNumBinEl++;
  }

  /**
   * Remove the element at the given index
   * 
   * @attention Constant time!
   * @param idx element index
   * @internal @return index of the binned element that was removed
   */
  UINTEGER PSSACR_Bin::remove_at(UINTEGER idx)
  {
    if( (idx + 1) == unNumBinEl)
      --unNumBinEl;
    else
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(idx < unNumBinEl)
#endif
      arunBinEl[idx] = arunBinEl[--unNumBinEl];
#ifndef PSSALIB_NO_BOUNDS_CHECKS
    else
      throw std::runtime_error("PSSACR_Bin::remove_at() - index out of range.");
#endif

    return arunBinEl[idx]; // element number to be updated
  }

  /**
   * Retrieve element at the given index
   * 
   * @param idx index of the element
   * @return binned propensity index
   */
  UINTEGER PSSACR_Bin::get_at(UINTEGER idx) const
  {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(idx >= unNumBinEl)
        throw std::runtime_error("PSSACR_Bin::get_at() - index out of range.");
      else
#endif
        return arunBinEl[idx];
  }

  /**
   * Returns number of entries in index list
   * @return number of entries in index list
   */
  UINTEGER PSSACR_Bin::size() const
  {
    return unNumBinEl;
  }

  ////////////////////////////////
  // Bin storage & mapping

  //! Default constructor
  PSSACR_Bins::PSSACR_Bins()
    : binVals (NULL)
    , unVals(0)
  {
    // Goggle dense_hash_map specific
#ifdef __USE_GOOGLE_HASH_MAP
    mapBins.set_empty_key(0);
#endif
  }

  //! Destructor
  PSSACR_Bins::~PSSACR_Bins()
  {
    clear();
  }

  /**
   * @brief Updates value in a bin
   *
   * @param bin_no_new new bin number for the given index. Zero denotes a deletion.
   * @param idx item index (e.g., row index in group sum array).
   * @param val new value to be set for item at index idx.
   */
  void PSSACR_Bins::updateValue(UINTEGER bin_no_new, UINTEGER idx, REAL val)
  {
    if(val > 0.0)
    {
      if(0 == binVals[idx].bin_no)
      {
        // insert
        it = mapBins.find(bin_no_new);
        if(mapBins.end() == it)
        {
          PSSACR_Bin bin;
          insRes = mapBins.insert(MAP_BINS::value_type(bin_no_new, bin));

          it = insRes.first;
          it->second.dBinSum = val;
          it->second.resize(unVals);

          binVals[idx].idx = it->second.push_back(idx);
          binVals[idx].bin_no = bin_no_new;
          binVals[idx].val = val;
        }
        else
        {
          it->second.dBinSum += val;

          binVals[idx].idx = it->second.push_back(idx);
          binVals[idx].bin_no = bin_no_new;
          binVals[idx].val = val;
        }
      }
      else if(bin_no_new == binVals[idx].bin_no)
      {
        // update
        it = mapBins.find(binVals[idx].bin_no);
        if(mapBins.end() == it)
        {
          throw std::runtime_error("PSSACR_Bin::updateValue() - bin does not exist");
          //std::cerr << "Bin does not exist!\nBin = " << binVals[idx].bin_no
          //  << " : idx = " << idx << " : val = " << val << std::endl;
        }
        else
        {
          it->second.dBinSum +=  val - binVals[idx].val;

          binVals[idx].val = val;
        }
      }
      else
      {
        // update & move
        it = mapBins.find(binVals[idx].bin_no);
        if(mapBins.end() == it)
        {
          throw std::runtime_error("PSSACR_Bin::updateValue() - bin does not exist");
          //std::cerr << "Bin does not exist!\nBin = " << binVals[idx].bin_no
          //  << " : idx = " << idx << " : val = " << val << std::endl;
        }
        else
        {
          // remove at old position & update the index of swapped element
          binVals[it->second.remove_at(binVals[idx].idx)].idx = binVals[idx].idx;

          it->second.dBinSum -= binVals[idx].val;

          it = mapBins.find(bin_no_new);
          if(mapBins.end() == it)
          {
            // insert new bin
            PSSACR_Bin bin;
            insRes = mapBins.insert(MAP_BINS::value_type(bin_no_new, bin));

            it = insRes.first;
            it->second.dBinSum = val;
            it->second.resize(unVals);

            binVals[idx].idx = it->second.push_back(idx);
          }
          else
          {
            it->second.dBinSum +=  val;

            binVals[idx].idx = it->second.push_back(idx);
          }

          binVals[idx].bin_no = bin_no_new;
          binVals[idx].val = val;
        }
      }
    }
    else
    {
      it = mapBins.find(binVals[idx].bin_no);
      if(mapBins.end() != it)
      {
        // delete empty bin entry & update the index of swapped element
        binVals[it->second.remove_at(binVals[idx].idx)].idx = binVals[idx].idx;
        it->second.dBinSum -= binVals[idx].val;

        binVals[idx].bin_no = 0;
        binVals[idx].val = 0.0;
        binVals[idx].idx = -1;
      }
    }
  }
}  } // close namespaces pssalib and datamodel
