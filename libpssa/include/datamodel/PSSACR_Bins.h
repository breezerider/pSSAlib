/**
 * @file PSSACR_Bins.h
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
 * Binning is created along with an instance of @c PSSACR_Bin class.
 * @c PSSACR_Bin class implements the @c updateValue function
 * that handles all the updates necessary when a value of
 * a binned quantity got changed. Check out the source for details.
 */

#ifndef PSSALIB_DATAMODEL_PSSACR_BINS_H_
#define PSSALIB_DATAMODEL_PSSACR_BINS_H_

#include "../stdheaders.h"
#include "../typedefs.h"

namespace pssalib
{
namespace datamodel
{
  ////////////////////////////////
  //! \class PSSACR_Bin
  //! \brief Bin class
  class PSSACR_Bin
  {
    ////////////////////////////////
    // Attributes
    public:
      //! Sum of bin elements
      REAL dBinSum;
      //! Indices of binned elements in the original array
      UINTEGER *arunBinEl;

      //! \internal auxiliary variables:
      UINTEGER unNumBinEl;
      //! \internal auxiliary variables:
      UINTEGER unCapBinEl;

    ////////////////////////////////
    // Constructors
    public:
      // Constructor
      PSSACR_Bin();

      // Copy constructor
      PSSACR_Bin(const PSSACR_Bin &b);

      // Destructor
      virtual ~PSSACR_Bin();

    ////////////////////////////////
    // Methods
    public:
      // Add an index to the list
      UINTEGER push_back(UINTEGER el);
      // Remove the element at the given index
      UINTEGER remove_at(UINTEGER idx);
      // Retrieve element at the given index
      UINTEGER get_at(UINTEGER idx) const;

      // Returns number of entries in index list
      UINTEGER size() const;
      // Resize the bin
      void resize(UINTEGER n);
      // Clear the bin
      void clear();
  };

  ////////////////////////////////
  //! \class PSSACR_Bins
  //! \brief Bins storage & mapping
  class PSSACR_Bins
  {
  ////////////////////////////////
  // Typedefs
  public:
    typedef HASHMAP<UINTEGER, PSSACR_Bin> MAP_BINS;
    typedef HASHMAP<UINTEGER, PSSACR_Bin>::iterator BINS_ITER;
    typedef HASHMAP<UINTEGER, PSSACR_Bin>::const_iterator CONST_BINS_ITER;
    typedef std::pair<BINS_ITER, BINS_ITER> PAIR_BINS_ITER;
    typedef std::pair<CONST_BINS_ITER, CONST_BINS_ITER> CONST_PAIR_BINS_ITER;

    //! Structure for fast access to bins
    typedef struct tagBinsVals
    {
      //! binned element value
      REAL val;
      //! mapping variable: Bin number
      UINTEGER bin_no;
      //! mapping variable: Index of the corresponding element in the bin
      UINTEGER idx;
      //! Default constructor
      tagBinsVals():
        val(0.0), bin_no(0), idx(0)
      {
        // do nothing
      };
    } BinVals;

  ////////////////////////////////
  // Attributes
  protected:
    MAP_BINS mapBins;
    BinVals *binVals;
    UINTEGER unVals;

    //! @internal Temporary variables
    BINS_ITER it;
    //! @internal Temporary variables
    std::pair<BINS_ITER, bool> insRes;

  ////////////////////////////////
  // Constructors
  public:
    // Constructor
    PSSACR_Bins();

    // Destructor
  virtual ~PSSACR_Bins();

    ////////////////////////////////
    // Methods
    public:

    /**
     * Get a pair of iterators for iterating across all bins
     * 
     * @param pit iterator pair
     */
  inline void getBins(CONST_PAIR_BINS_ITER &pit) const
    {
      pit.first = mapBins.begin();
      pit.second = mapBins.end();
    }

    //! Clear bins
  inline void clear()
    {
      if(mapBins.size() > 0)
        mapBins.clear();
      if(binVals != NULL)
      {
        delete [] binVals;
        binVals = NULL;
        unVals = 0;
      }
    }

    /**
     * Resizes the bins to ensuring they have
     * enough capacity to store N elements
     * 
     * @param N new bin size
     */
    inline void resize(UINTEGER N)
    {
      clear();
      mapBins.rehash(N);
      binVals = new BinVals[N];
      unVals = N;
    }

    // Updates value in a bin
    void updateValue(UINTEGER bin_no_new, UINTEGER idx, REAL val);

    //! Get value of bin with index <i>idx</i>
   inline REAL getValue(UINTEGER idx) const
    {
      if(idx < unVals)
        return binVals[idx].val;
      else
        throw std::runtime_error("PSSACR_Bins::getValue() - invalid index in arguments.");
    }
  };

}  } // close namespaces pssalib and datamodel

#endif /* PSSALIB_DATAMODEL_PSSACR_BINS_H_ */
