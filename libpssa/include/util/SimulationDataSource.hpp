/**
 * @file SimulationDataSource.hpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Feb 2017
 * @license The GNU LGPL v3 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Auxiliary classes for parsing pSSAlib data sets
 */

#include "../typedefs.h"
#include "MPIWrapper.h"

#ifndef PSSALIB_UTIL_SIMULATION_DATA_SOURCE_HPP_
#define PSSALIB_UTIL_SIMULATION_DATA_SOURCE_HPP_

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/**
 * @class OutputFormatter
 * @brief Interface for an output formatter
 */
class OutputFormatter
{
///////////////
// Data Types
public:
  //! Data dimensions
  enum DataDimensions
  {
    ddTime,
    ddSubvolume,
    ddSpecies
  };

///////////////
// Methods
public:

  /**
   * Map data dimension to a loop index
   * 
   * @param idx A dimensions of the dataset (can be time, subvolume or species)
   * @return corresponding loop index
   */
virtual BYTE mapIndex(DataDimensions dim) const = 0;

  /**
   * Get the data delimiter
   * 
   * @return character separator
   */
virtual const STRING::value_type * getDataSeparator() const = 0;

  /**
   * Get the file extension
   * 
   * @return character string representing filee extension
   */
virtual const STRING::value_type * getFileExtension() const = 0;

  /**
   * Check if splitting accross multiple files is desired
   * 
   * @return @c true if data is to be split along the first ordinal index, @c false otherwise
   */
virtual bool isSplittable() const = 0;

  /**
   * Format header for a given ordinal index
   * 
   * @param idx Ordinal index of the loop
   * @return character string representing a header
   */
virtual const STRING::value_type * getHeader(BYTE idx, UINTEGER pos) const = 0;

  /**
   * Format footer for a given ordinal index
   * 
   * @param idx Ordinal index of the loop
   * @return character string representing a footer
   */
virtual const STRING::value_type * getFooter(BYTE idx, UINTEGER pos) const = 0;
};

/**
 * @class CSVOutputFormatter
 * @brief A CSV output formatter
 */
class CSVOutputFormatter : public OutputFormatter
{
///////////////
// Attributes
protected:
  //!
  STRING m_strHeader;

  REAL   m_dTimeBegin,
         m_dTimeStep;

///////////////
// Constructors
public:
  //! Constructor
  CSVOutputFormatter(const STRING & header = STRING(),
                     const REAL dt = 0.0,
                     const REAL tb = 0.0)
    : m_strHeader(header)
    , m_dTimeBegin(tb)
    , m_dTimeStep(dt)
  {
    // Do nothing
  }

///////////////
// Methods
public:

  /**
   * @copydoc OutputFormatter::mapIndex(DataDimensions)
   */
virtual BYTE mapIndex(DataDimensions dim) const
  {
    switch(dim)
    {
    case ddTime:
      return 1;
    case ddSubvolume:
      return 0;
    case ddSpecies:
      return 2;
    default:
      return BYTE(-1);
    }
  }

  /**
   * @copydoc OutputFormatter::getDataSeparator()
   */
virtual const STRING::value_type * getDataSeparator() const
  {
    return ",";
  }

  /**
   * @copydoc OutputFormatter::getFileExtension()
   */
virtual const STRING::value_type * getFileExtension() const
  {
    return "csv";
  }

  /**
   * @copydoc OutputFormatter::isSplittable()
   */
virtual bool isSplittable() const
  {
    return true;
  }

  /**
   * @copydoc OutputFormatter::getHeader(BYTE)
   */
virtual const STRING::value_type * getHeader(BYTE idx, UINTEGER pos) const
  {
    if(0 == idx)
      return m_strHeader.c_str();
    else if((1 == idx)&&(0.0 < m_dTimeStep))
    {
      static STRING        strResult;
      static STRINGSTREAM  ssTemp;

      ssTemp << (m_dTimeBegin + pos*m_dTimeStep) << getDataSeparator();
      strResult = ssTemp.str();
      ssTemp.str(STRING());
      return strResult.c_str();
    }
    else
      return NULL;
  }

  /**
   * @copydoc OutputFormatter::getFooter(BYTE)
   */
virtual const STRING::value_type * getFooter(BYTE idx, UINTEGER pos) const
  {
    if(2 == idx)
      return "\n";
    else
      return NULL;
  }
};

/**
 * @class VTKOutputFormatter
 * @brief A VTK output formatter
 */
class VTKOutputFormatter : public OutputFormatter
{
///////////////
// Attributes
protected:
  //! 
  BYTE                m_uDims;
  const UINTEGER *    m_arDims;

  //!
  std::vector<STRING> m_arSpeciesIds;

///////////////
// Constructors
public:
  //! Constructor
  VTKOutputFormatter(BYTE uDims, const UINTEGER * pDims, const std::vector<STRING> & arIds)
    : m_uDims(uDims)
    , m_arDims(pDims)
    , m_arSpeciesIds(arIds)
  {
    if(NULL == m_arDims)
      m_uDims = 0;
  }

///////////////
// Methods
public:
  /**
   * @copydoc OutputFormatter::mapIndex(DataDimensions)
   */
virtual BYTE mapIndex(DataDimensions dim) const
  {
    switch(dim)
    {
    case ddTime:
      return 0;
    case ddSubvolume:
      return 2;
    case ddSpecies:
      return 1;
    default:
      return BYTE(-1);
    }
  }

  /**
   * @copydoc OutputFormatter::getDataSeparator()
   */
virtual const STRING::value_type * getDataSeparator() const
  {
    return "\n";
  }

  /**
   * @copydoc OutputFormatter::getFileExtension()
   */
virtual const STRING::value_type * getFileExtension() const
  {
    return "vtk";
  }

  /**
   * @copydoc OutputFormatter::isSplittable()
   */
virtual bool isSplittable() const
  {
    return true;
  }

  /**
   * @copydoc OutputFormatter::getHeader(BYTE)
   */
virtual const STRING::value_type * getHeader(BYTE idx, UINTEGER pos) const
  {
    static STRING strResult;
    STRINGSTREAM ssTemp;
    if(0 == idx)
    {
      UINTEGER subvolumes = 1;
      for(BYTE di = 0; di < m_uDims; ++di)
        subvolumes *= m_arDims[di];

      ssTemp << "# vtk DataFile Version 3.0\n";
      ssTemp << "Partial propensity vtk frame: " << pos << " volumes: " << subvolumes << "\n";

      ssTemp << "ASCII\n";
      ssTemp << "DATASET STRUCTURED_POINTS\n";

      ssTemp << "DIMENSIONS " << UINTEGER(m_uDims > 0 ? m_arDims[0]+1 : 1) << " " << UINTEGER(m_uDims > 1 ? m_arDims[1]+1 : 1) << " " << UINTEGER(m_uDims > 2 ? m_arDims[2]+1 : 1) << "\n";

      ssTemp << "ORIGIN 0 0 0\n";
      ssTemp << "SPACING 1 1 1\n";
      ssTemp << "CELL_DATA " << subvolumes << "\n";

      strResult = ssTemp.str();
      return strResult.c_str();
    }
    else if(1 == idx)
    {
      if(m_arSpeciesIds.size() > pos)
        ssTemp << "SCALARS " << m_arSpeciesIds[pos];
      else
        ssTemp << "SCALARS species" << pos;
      ssTemp << " unsigned_int 1\n";
      ssTemp << "LOOKUP_TABLE default\n";

      strResult = ssTemp.str();
      return strResult.c_str();
    }

    return NULL;
  }

  /**
   * @copydoc OutputFormatter::getFooter(BYTE)
   */
virtual const STRING::value_type * getFooter(BYTE idx, UINTEGER pos) const
  {
    if(2 == idx)
      return "\n";
    else if(1 == idx)
      return "\r";
    else
      return NULL;
  }
};

/**
 * @class GnuplotOutputFormatter
 * @brief A GNUPlot output formatter
 */
class GnuplotOutputFormatter : public OutputFormatter
{
///////////////
// Attributes
protected:
  //!
  STRING m_strHeader;

  REAL   m_dTimeBegin,
         m_dTimeStep;

///////////////
// Constructors
public:
  //! Constructor
  GnuplotOutputFormatter(const STRING & header = STRING(),
                     const REAL dt = 0.0,
                     const REAL tb = 0.0)
    : m_strHeader(header)
    , m_dTimeBegin(tb)
    , m_dTimeStep(dt)
  {
    // Do nothing
  }

///////////////
// Methods
public:

  /**
   * @copydoc OutputFormatter::mapIndex(DataDimensions)
   */
virtual BYTE mapIndex(DataDimensions dim) const
  {
    switch(dim)
    {
    case ddTime:
      return 1;
    case ddSubvolume:
      return 0;
    case ddSpecies:
      return 2;
    default:
      return BYTE(-1);
    }
  }

  /**
   * @copydoc OutputFormatter::getDataSeparator()
   */
virtual const STRING::value_type * getDataSeparator() const
  {
    return " ";
  }

  /**
   * @copydoc OutputFormatter::getFileExtension()
   */
virtual const STRING::value_type * getFileExtension() const
  {
    return "gnuplot";
  }

  /**
   * @copydoc OutputFormatter::isSplittable()
   */
virtual bool isSplittable() const
  {
    return false;
  }

  /**
   * @copydoc OutputFormatter::getHeader(BYTE)
   */
virtual const STRING::value_type * getHeader(BYTE idx, UINTEGER pos) const
  {
    if(0 == idx)
      return m_strHeader.c_str();
    else if((1 == idx)&&(0.0 < m_dTimeStep))
    {
      static STRING        strResult;
      static STRINGSTREAM  ssTemp;

      ssTemp << (m_dTimeBegin + pos*m_dTimeStep) << getDataSeparator();
      strResult = ssTemp.str();
      ssTemp.str(STRING());
      return strResult.c_str();
    }
    else
      return NULL;
  }

  /**
   * @copydoc OutputFormatter::getFooter(BYTE)
   */
virtual const STRING::value_type * getFooter(BYTE idx, UINTEGER pos) const
  {
    if(2 == idx)
      return "\n";
    if(0 == idx)
      return "\n \n \n";
    else
      return NULL;
  }
};

/**
 * @class MultIterator
 * @brief An aggregated multi-dimensional iterator
 */
class MultIterator
{
private:
  const UINTEGER unDoNotIterate;
  
protected:
  UINTEGER **indexes;
  
  BYTE idxLoop, numLoops;

public:
  MultIterator(BYTE n)
   : unDoNotIterate(0)
   , indexes(NULL)
   , idxLoop(0)
   , numLoops(n)
  {
    indexes = new UINTEGER*[numLoops*4];
    memset(indexes, (char)0, sizeof(UINTEGER)*numLoops*4);
  }

  ~MultIterator()
  {
    delete [] indexes;
  }

public:
inline BYTE idx()
  {
    return idxLoop;
  }

inline BYTE len()
  {
    const void * ptrDoNotIterate = &unDoNotIterate;
    if(ptrDoNotIterate == indexes[idxLoop + numLoops])
      return (*indexes[idxLoop + 3*numLoops]);
    else
      return indexes[idxLoop + 3*numLoops] - indexes[idxLoop + 2*numLoops];
  }

inline void advance()
  {
    ++idxLoop; reset();
  }

inline void retreat()
  {
    --idxLoop;
  }

  UINTEGER pos()
  {
    return (*indexes[idxLoop]);
  }

inline void init()
  {
    for(idxLoop = 0; idxLoop < numLoops; ++idxLoop)
      reset();
    idxLoop = 0;
  }

  bool set(BYTE idxLoop, UINTEGER * pOut, const UINTEGER * pMin, const UINTEGER * pMax, bool iterate = true)
  {
    const void * ptrDoNotIterate = &unDoNotIterate;
    if(idxLoop >= numLoops)
      return false;
    indexes[idxLoop] = pOut; indexes[idxLoop + 2*numLoops] = (UINTEGER *)pMin; indexes[idxLoop + 3*numLoops] = (UINTEGER *)pMax;
    if(!iterate)
      indexes[idxLoop + numLoops] = (UINTEGER*)ptrDoNotIterate;
    return true;
  }

  void increment()
  {
    const void * ptrDoNotIterate = &unDoNotIterate;
    if(ptrDoNotIterate == indexes[idxLoop + numLoops])
      ++(*indexes[idxLoop]);
    else
    {
      ++indexes[idxLoop + numLoops];
      *indexes[idxLoop] = *indexes[idxLoop + numLoops];
    }
  }

  bool good()
  {
    const void * ptrDoNotIterate = &unDoNotIterate;
    if(ptrDoNotIterate == indexes[idxLoop + numLoops])
      return (*indexes[idxLoop]) < (*indexes[idxLoop + 3*numLoops]);
    else
      return indexes[idxLoop + numLoops] != indexes[idxLoop + 3*numLoops];
  }

  void reset()
  {
    const void * ptrDoNotIterate = &unDoNotIterate;
    if(ptrDoNotIterate == indexes[idxLoop + numLoops])
      *indexes[idxLoop] = *indexes[idxLoop + 2*numLoops];
    else
    {
      indexes[idxLoop + numLoops] = indexes[idxLoop + 2*numLoops];
      *indexes[idxLoop] = *indexes[idxLoop + numLoops];
    }
  }
};

/**
 * @class SimulationDataSource
 * @brief Auxiliary class for parsing a pSSAlib data set
 */
class SimulationDataSource
{
///////////////
// Attributes
protected:
  // Data source specs
  UINTEGER m_unTimePoints, //!< Number of time points
           m_unSpecies,    //!< Number of species
           m_unSubvolumes; //!< Number of subvolumes

  //! Data
  REAL *m_arData;

///////////////
// Constructors
public:

  SimulationDataSource(UINTEGER rows = 0,
                       UINTEGER cols = 0,
                       UINTEGER vols = 0
                      )
    : m_unTimePoints(rows)
    , m_unSpecies(cols)
    , m_unSubvolumes(vols)
    , m_arData(NULL)
  {
    if(m_unTimePoints*m_unSpecies > 0)
    {
      if(0 == m_unSubvolumes) m_unSubvolumes = 1;

      alloc();
    }
  };

  ~SimulationDataSource()
  {
    free();
  };

///////////////
// Methods
protected:

  /**
   * Allocate memory storage for data
   * 
   * @return @true on success, @false otherwise
   */
inline bool alloc()
  {
    if(m_unTimePoints*m_unSpecies*m_unSubvolumes > 0)
    {
      try
      {
        m_arData = new REAL[m_unTimePoints*m_unSpecies*m_unSubvolumes];
      }
      catch(std::bad_alloc &e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : could not allocate data storage, reason: '" << e.what() << "'.\n";
        return false;
      }

      memset(m_arData, (char)0, m_unSpecies*sizeof(void*));

      return true;
    }
    return false;
  }

  /**
   * Release the memory storage
   */
inline void free()
  {
    if(NULL != m_arData)
    {
      delete [] m_arData;
      m_unTimePoints = m_unSpecies = m_unSubvolumes = 0;
    }
    m_arData = NULL;
  }

///////////////
// Methods
public:

  /**
   * Clean up
   */
  inline void clear()
  {
    free();
  }

  /**
   * Get number of time points
   */
inline UINTEGER getTimePoints() const
  {
    return m_unTimePoints;
  }

  /**
   * Get number of species
   */
inline UINTEGER getSpecies() const
  {
    return m_unSpecies;
  }

  /**
   * Get number of subvolumes
   */
inline UINTEGER getSubvolumes() const
  {
    return m_unSubvolumes;
  }

  /**
   * Get data entry at a given index
   * 
   * @param time Time index
   * @param species Species index
   * @param subvol Subvolume index
   * @return Data entry at given position
   */
  REAL & at(UINTEGER time, UINTEGER species, UINTEGER subvol = 0)
  {
    static REAL empty;

    if((time >= m_unTimePoints)||(species >= m_unSpecies)||(subvol >= m_unSubvolumes))
    {
      PSSALIB_MPI_CERR_OR_NULL << "index out of bounds (" << time << ", "
        << species << ", " << subvol << ")\n";
      empty = REAL();
      return empty;
    }
    return m_arData[time*m_unSpecies*m_unSubvolumes + subvol*m_unSpecies + species];
  }

  /**
   * Load a simulation data set from a file
   * 
   * @param filePath Path to data set file
   * @param rangeTime Initial and final time points to be processed
   * @param rangeSubvolumes Pointer to an array containing subvolume indexes
   * @param rangeSpecies Pointer to an array containing species indexes
   * @return @c true on success, @c false otherwise
   */
virtual bool load(const STRING & filePath,
                  const std::pair< UINTEGER, UINTEGER > & rangeTime = std::pair< UINTEGER, UINTEGER >(),
                  const std::pair< const UINTEGER *, const UINTEGER * > & rangeSpecies = std::pair< const UINTEGER *, const UINTEGER * >(),
                  const std::pair< const UINTEGER *, const UINTEGER * > & rangeSubvolumes = std::pair< const UINTEGER *, const UINTEGER * >()
                 )
{
  // Argument checks
  if(filePath.empty())
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error : an empty file path provided.\n";
    return false;
  }

  bool result = false;
  FILESTREAMBUFFER fsbData;
  if(!fsbData.open(filePath.c_str(), std::ios_base::in))
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error : could not open file '" << filePath << "'\n";
  }
  else
  {
    ISTREAM isData(&fsbData);
    result = load(isData, rangeTime, rangeSpecies, rangeSubvolumes);
  }

  fsbData.close();

  return result;
}

  /**
   * Load a simulation data set from a stream
   * 
   * @param filePath Path to data set file
   * @param rangeTime Initial and final time points to be processed
   * @param rangeSubvolumes Pointer to an array containing subvolume indexes
   * @param rangeSpecies Pointer to an array containing species indexes
   * @return @c true on success, @c false otherwise
   */
virtual bool load(ISTREAM & isData,
                  const std::pair< UINTEGER, UINTEGER > & rangeTime = std::pair< UINTEGER, UINTEGER >(),
                  const std::pair< const UINTEGER *, const UINTEGER * > & rangeSpecies = std::pair< const UINTEGER *, const UINTEGER * >(),
                  const std::pair< const UINTEGER *, const UINTEGER * > & rangeSubvolumes = std::pair< const UINTEGER *, const UINTEGER * >()
                 )
  {
    // Argument checks
    if(!isData.good())
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : input stream is invalid\n";
      return false;
    }
    
    if(rangeTime.second < rangeTime.first)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid temporal range (" << rangeTime.first << ", " << rangeTime.second << ")\n";
      return false;
    }

    if(((NULL != rangeSpecies.first)&&(NULL == rangeSpecies.second))||
       ((NULL == rangeSpecies.first)&&(NULL != rangeSpecies.second))||
       ((NULL != rangeSpecies.first)&&(NULL != rangeSpecies.second)&&
        (rangeSpecies.first == rangeSpecies.second))||
       (rangeSpecies.first > rangeSpecies.second))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid species range\n";
      return false;
    }
    
    if(((NULL != rangeSubvolumes.first)&&(NULL == rangeSubvolumes.second))||
       ((NULL == rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second))||
       ((NULL != rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second)&&
        (rangeSubvolumes.first == rangeSubvolumes.second))||
       (rangeSubvolumes.first > rangeSubvolumes.second))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid subvolume range\n";
      return false;
    }

    // Time
    if(rangeTime.second > 0)
    {
      m_unTimePoints = rangeTime.second - rangeTime.first;
    }

    // Subvolumes
    std::set<UINTEGER> setIdxSubvolumes;
    if((NULL != rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second))
    {
      setIdxSubvolumes.insert(rangeSubvolumes.first, rangeSubvolumes.second);
      m_unSubvolumes = setIdxSubvolumes.size();
    }

    // Species
    std::set<UINTEGER> setIdxSpecies;
    if((NULL != rangeSpecies.first)&&(NULL != rangeSpecies.second))
    {
      setIdxSpecies.insert(rangeSpecies.first, rangeSpecies.second);
      m_unSpecies = setIdxSpecies.size();
    }

    if(0 == m_unTimePoints)
    {
      // seek towards the beginning
      isData.seekg(0, isData.beg);

      // new lines will be skipped unless we stop it from happening
      isData.unsetf(std::ios_base::skipws);

      // count the newlines with an algorithm specialized for counting
      m_unTimePoints = std::count(
          std::istream_iterator<STRING::value_type>(isData),
          std::istream_iterator<STRING::value_type>(), 
          '\n');
      
      if(0 == m_unTimePoints)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : empty data set\n";
        return false;
      }
      
      // restore the flag & clear the bits
      isData.setf(std::ios_base::skipws);
      isData.clear();

      // check if file ends with a line break and decrease line number accordingly
      isData.seekg(-1, isData.end);
      if(isData.peek() != '\n')
        ++m_unTimePoints;
      
      // seek towards the beginning
      isData.seekg(0, isData.beg);
    }

    UINTEGER unCurrTimePoint = 0, unCurrSubvolume = 0, unCurrSpecies = 0;
    UINTEGER timePoint = 0, subVol = 0, species = 0,
             unSpecies = 0, unSubvolumes = 0;

    try
    {
      static boost::char_separator<STRING::value_type> 
        sep_subvol(PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER), sep_col(PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER);
      STRING strTemp;

      for(;unCurrTimePoint < m_unTimePoints; ++timePoint)
      {
        if(!isData.good()||!getline(isData, strTemp))
          throw std::runtime_error("unexpected end of file");

        boost::algorithm::trim(strTemp);
        if(strTemp.empty())
          throw std::runtime_error("file cannot contain blank lines");

        if(timePoint < rangeTime.first)
        {
          continue;
        }

        unSubvolumes = 0;
        for(STRING::size_type pos = STRING::size_type(0); ; ++unSubvolumes, ++pos)
        {
          pos = strTemp.find(PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER, pos);
          if(STRING::npos == pos)
            break;
        }
        if((strTemp.length()-1) != strTemp.rfind(PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER))
          ++unSubvolumes;
        if(0 == unSubvolumes) unSubvolumes = 1;
        if(0 == m_unSubvolumes)
          m_unSubvolumes = unSubvolumes;
        else if(unSubvolumes < m_unSubvolumes)
          throw std::runtime_error("not enought subvolumes");

        boost::tokenizer<
          boost::char_separator<STRING::value_type> >
          tok_subvol(strTemp, sep_subvol);

        boost::tokenizer<
          boost::char_separator<STRING::value_type> >::iterator
          itTokSubVol = tok_subvol.begin();

        std::set<UINTEGER>::iterator idxSubvolumesIt = setIdxSubvolumes.begin(); unCurrSubvolume = 0;
        for(subVol = 0; subVol < unSubvolumes; ++subVol)
        {
          if(0 < setIdxSubvolumes.size())
          {
            if(setIdxSubvolumes.end() != idxSubvolumesIt)
            {
              if(subVol < *idxSubvolumesIt)
              {
                if(tok_subvol.end() == itTokSubVol)
                  throw std::runtime_error("inconsistent number of subvolume");

                itTokSubVol++;
                continue;
              }
              idxSubvolumesIt++;
            }
            else
              break; // done
          }

          unSpecies = 0;
          for(STRING::size_type pos = STRING::size_type(0); ; ++unSpecies, ++pos)
          {
            pos = (*itTokSubVol).find(PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER, pos);
            if(STRING::npos == pos)
              break;
          }
          if(((*itTokSubVol).length()-1) != (*itTokSubVol).rfind(PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER))
            ++unSpecies;
          if(0 == unSpecies) unSpecies = 1;
          if(0 == m_unSpecies)
            m_unSpecies = unSpecies;
          else if(unSpecies < m_unSpecies)
            throw std::runtime_error("not enought species");

          boost::tokenizer<
            boost::char_separator<STRING::value_type> >
            tok_cols(*itTokSubVol, sep_col);

          boost::tokenizer<
            boost::char_separator<STRING::value_type> >::iterator
            itTokCols = tok_cols.begin();

          std::set<UINTEGER>::iterator idxSpeciesIt = setIdxSpecies.begin(); unCurrSpecies = 0;
          for(species = 0; (species < unSpecies)&&(unCurrSpecies < m_unSpecies); ++species)
          {
            if(0 < setIdxSpecies.size())
            {
              if(setIdxSpecies.end() != idxSpeciesIt)
              {
                if(species < *idxSpeciesIt)
                {
                  if(tok_cols.end() == itTokCols)
                    throw std::runtime_error("inconsistent species number accross subvolumes");

                  itTokCols++;
                  continue;
                }
                idxSpeciesIt++;
              }
              else
                break; // done
            }
            
            if((NULL == m_arData)&&!alloc())
              throw std::runtime_error("could not allocate memory");

            at(unCurrTimePoint, unCurrSpecies, unCurrSubvolume) = boost::lexical_cast<REAL>(*itTokCols);
            
            ++unCurrSpecies; itTokCols++;
          }
          
          ++unCurrSubvolume; itTokSubVol++;
        }
        ++unCurrTimePoint;
      }
    }
    catch(std::runtime_error & e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : '" << e.what()
        << "' when processing input stream on line "
        << timePoint+1 << " at subvolume " << subVol+1
        << " species " << species+1 << ".\n";
      return false;
    }
    catch(boost::bad_lexical_cast & e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : '" << e.what()
        << "' when processing input stream on line "
        << timePoint+1 << " at subvolume " << subVol+1 
        << " species " << species+1 << ".\n";
      return false;
    }
    catch(...)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : unexpected "
        "exception when processing input stream on line "
        << timePoint+1 << " at subvolume " << subVol+1 
        << " species " << species+1 << ".\n";
      return false;
    }

    return true;
  }

  /**
   * Store a simulation data setusing a given format
   * 
   * @param filePath Path to output destination
   * @param fmt An @c OutputFormatter object to format the output
   * @param rangeTime Initial and final time points to be processed
   * @param rangeSubvolumes Pointer to an array containing subvolume indexes
   * @param rangeSpecies Pointer to an array containing species indexes
   * @return @c true on success, @c false otherwise
   */
virtual bool store(const STRING & filePath, const OutputFormatter & fmt = CSVOutputFormatter(),
                   const std::pair< UINTEGER, UINTEGER > & rangeTime = std::pair< UINTEGER, UINTEGER >(),
                   const std::pair< const UINTEGER *, const UINTEGER * > & rangeSpecies = std::pair< const UINTEGER *, const UINTEGER * >(),
                   const std::pair< const UINTEGER *, const UINTEGER * > & rangeSubvolumes = std::pair< const UINTEGER *, const UINTEGER * >()
                 )
  {
    // Argument checks
    if(filePath.empty())
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : an empty file path provided.\n";
      return false;
    }

    if((rangeTime.second < rangeTime.first)||(rangeTime.second >= m_unTimePoints))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid temporal range (" << rangeTime.first << ", " << rangeTime.second << ") with only " << m_unTimePoints << " time points available.\n";
      return false;
    }

    if(((NULL != rangeSpecies.first)&&(NULL == rangeSpecies.second))||
       ((NULL == rangeSpecies.first)&&(NULL != rangeSpecies.second))||
       ((NULL != rangeSpecies.first)&&(NULL != rangeSpecies.second)&&
        (rangeSpecies.first == rangeSpecies.second))||
       (rangeSpecies.first > rangeSpecies.second))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid species range\n";
      return false;
    }

    if(((NULL != rangeSubvolumes.first)&&(NULL == rangeSubvolumes.second))||
       ((NULL == rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second))||
       ((NULL != rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second)&&
        (rangeSubvolumes.first == rangeSubvolumes.second))||
       (rangeSubvolumes.first > rangeSubvolumes.second))
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : invalid subvolume range\n";
      return false;
    }

    // Vars
    const BYTE numLoops = 3;
    UINTEGER timePoint = 0, subVol = 0, species = 0;
    const UINTEGER zero = 0;
    STRING strFilePathPattern(filePath), strFilePath;
    
    try
    {
      // Iterator
      MultIterator its(numLoops);

      // Time
      UINTEGER unEndTimePoint = ((rangeTime.second > 0) && (rangeTime.second < m_unTimePoints)) ? rangeTime.second : m_unTimePoints;
      its.set(fmt.mapIndex(OutputFormatter::ddTime), &timePoint, &rangeTime.first, &unEndTimePoint, false);

      // Subvolumes
      boost::scoped_array<UINTEGER> arSubvolumes;
      if((NULL != rangeSubvolumes.first)&&(NULL != rangeSubvolumes.second))
      {
        std::size_t len = rangeSubvolumes.second - rangeSubvolumes.first;
        arSubvolumes.reset(new UINTEGER[len]);
        UINTEGER * end = std::copy(rangeSubvolumes.first, rangeSubvolumes.second, arSubvolumes.get());
        //std::sort(arSubvolumes.get(), end);
        while((*end) >= m_unSubvolumes) --end;
        
        its.set(fmt.mapIndex(OutputFormatter::ddSubvolume), &subVol, arSubvolumes.get(), end, true);
      }
      else
        its.set(fmt.mapIndex(OutputFormatter::ddSubvolume), &subVol, &zero, &m_unSubvolumes, false);

      // Species
      boost::scoped_array<UINTEGER> arSpecies;
      if((NULL != rangeSpecies.first)&&(NULL != rangeSpecies.second))
      {
        std::size_t len = rangeSpecies.second - rangeSpecies.first;
        arSpecies.reset(new UINTEGER[len]);
        UINTEGER * end = std::copy(rangeSpecies.first, rangeSpecies.second, arSpecies.get());
        //std::sort(arSpecies.get(), end);
        while((*end) >= m_unSpecies) --end;
        
        its.set(fmt.mapIndex(OutputFormatter::ddSpecies), &species, arSpecies.get(), end, true);
      }
      else
        its.set(fmt.mapIndex(OutputFormatter::ddSpecies), &species, &zero, &m_unSpecies, false);
  
      // initialize
      its.init();

      // I/O
      FILESTREAMBUFFER fsbData;
      OSTREAM osData(&fsbData);

      // check if extension is part of the file path
      std::size_t lenExt = strlen(fmt.getFileExtension()) + 1;
      bool bExtPresent = (strFilePathPattern.length() > lenExt ?
        (0 == strFilePathPattern.compare(strFilePathPattern.length() - lenExt,
                                         lenExt, STRING(".") + fmt.getFileExtension()))
        : false);

      // check if a placeholder for sequence index is present
      bool bSeqPresent = (STRING::npos != strFilePathPattern.find("%i")); 

      if(fmt.isSplittable()&&(0 < its.len()))
      {
        if(!bSeqPresent)
        {
          // ... and if not present add it
          if(!bExtPresent)
            strFilePathPattern.append("_%i.").append(fmt.getFileExtension());
          else
            strFilePathPattern.insert(strFilePathPattern.length() - lenExt, "_%i");
        }

        // test path
        strFilePath = (BOOSTFORMAT(strFilePathPattern) % 0).str();
        if(NULL == fsbData.open(strFilePath.c_str(), std::ios_base::out))
          return false;
        fsbData.close();
      }
      else
      {
        if(bSeqPresent)
          strFilePath = (BOOSTFORMAT(strFilePathPattern) % 0).str();
        else
          strFilePath = strFilePathPattern;

        if(!bExtPresent)
          strFilePath.append(".").append(fmt.getFileExtension());
        if(NULL == fsbData.open(strFilePath.c_str(), std::ios_base::out))
          return false;
      }

      // Main loop
      while(true)
      {
        // check termination critereon
        if(its.good())
        {
          // carry on
          if(fmt.isSplittable()&&(0 == its.idx())&&(0 < its.len()))
          {
            fsbData.close();
            strFilePath = (BOOSTFORMAT(strFilePathPattern) % its.pos()).str();
            if(NULL == fsbData.open(strFilePath.c_str(), std::ios_base::out))
              return false;
          }
          if(its.idx() != (numLoops - 1))
          {
            const STRING::value_type * p = fmt.getHeader(its.idx(), its.pos());
            if(NULL != p) osData << p;

            its.advance(); its.reset();
            if(its.idx() == (numLoops - 1))
            {
              const STRING::value_type * p = fmt.getHeader(its.idx(), its.pos());
              if(NULL != p) osData << p;
            }
            continue;
          }
          else
          {
            // output
            osData << m_arData[timePoint*m_unSubvolumes*m_unSpecies + subVol*m_unSpecies + species] << fmt.getDataSeparator();
          }
        }
        else
        {
          const STRING::value_type * p = fmt.getFooter(its.idx(), its.pos());
          if(NULL != p) osData << p;

          if(0 == its.idx())
            break;
          its.retreat();
        }
        its.increment();
      }
      fsbData.close();
    }
    catch(...)
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error : "
        << " when processing file '" << strFilePath 
        << "' on line " << timePoint+1 << " at subvolume " 
        << subVol+1 << " species " << species+1 << ".\n";
      return false;
    }

    return true;
  }
};

#endif /* PSSALIB_UTIL_SIMULATION_DATA_SOURCE_HPP_ */
