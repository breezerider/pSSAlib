/**
 * @file Subvolume_SPDM.hpp
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
 * Declares a container for subvolume variables required by the
 * Sorting Partial Propensity Direct Method (Ramaswamy, 2009)
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_SPDM_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_SPDM_HPP_

#include "../../stdheaders.h"
#include "JaggedMatrix.hpp"
#include "Subvolume_PDM.hpp"

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel_SPDM;

namespace detail
{
  /**
   * @copydoc Subvolume
   */
  class Subvolume_SPDM : public Subvolume_PDM
  {
  ////////////////////////////////
  // Friends
  public:
    friend class DataModel_SPDM;

  ////////////////////////////////
  // Data structures
  private:
    //! @internal sequence generator
    struct tagGenerateSequence
    {
      std::size_t index;

      tagGenerateSequence(std::size_t init = 0)
        : index(init)
      {
        // Do nothing
      };

      UINTEGER operator()()
      {
        return index++;
      };
    };

  ////////////////////////////////
  // Attributes
  protected:
    // Reactions
    //

    // Indexer variables
    std::size_t                 *m_IndexerRows;
    JaggedMatrix< std::size_t > m_IndexerCols;

  ////////////////////////////////
  // Constructors
  public:
    //! Constructor
    Subvolume_SPDM()
      : m_IndexerRows(NULL)
    {
      // Do nothing
    }

    //! Destructor
  virtual ~Subvolume_SPDM()
    {
      // Clean-up
      free_SPDM();
    }

  ////////////////////////////////
  // Methods
  private:

    /*
     * Free memory
     */
    void free_SPDM()
    {
      if(NULL != m_IndexerRows)
      {
        delete [] m_IndexerRows;
        m_IndexerRows = NULL;
      }
      m_IndexerCols.free();
    };

  ////////////////////////////////
  // Methods
  protected:
    /**
     * Reset all properties' values.
     */
  virtual void free()
    {
      // free memory
      free_SPDM();

      // call base class method
      Subvolume_PDM::free();
    };

    /**
     * @copydoc Subvolume::allocate(UINTEGER,UINTEGER,BYTE)
     */
  virtual void allocate(UINTEGER reactions, UINTEGER species, BYTE dims)
    {
      // call base class method
      Subvolume_PDM::allocate(reactions, species, dims);

      // allocate memory
      UINTEGER total_species = species + 1; // account for reservoir species
      m_IndexerRows = new std::size_t[total_species];
      memset(m_IndexerRows, (unsigned char)0, sizeof(std::size_t)*(total_species));
      m_IndexerCols.reserve(total_species, std::max(reactions / species, (UINTEGER)1));
    };

    /**
     * @copydoc Subvolume::clear(UINTEGER,UINTEGER)
     */
  virtual void clear(UINTEGER reactions, UINTEGER species)
    {
      memset(m_IndexerRows, (unsigned char)0, sizeof(std::size_t)*(species + 1)); // account for reservoir species
      m_IndexerCols.clear();

      // call base class method
      Subvolume_PDM::clear(reactions, species);
    };


  ////////////////////////////////
  // Methods
  public:

    /**
     * Reset mapping to original indexing 
     */
  inline void resetIndexing()
    {
      std::generate_n(m_IndexerRows, arPi.get_rows(), tagGenerateSequence());
      m_IndexerCols.resize(arPi.get_rows(), arPi.get_cols());
      for(std::size_t row = 0; row < arPi.get_rows(); ++row)
        for(std::size_t col = 0; col < arPi.get_cols(row); ++col)
          m_IndexerCols(row, col) = col;
    };

    /**
     * Map unsorted row index to the sorted index
     * 
     * @param index Unsorted row index
     * @return Mapped row index
     */
  inline std::size_t mapRowIndex(std::size_t i)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(i > unSpecies) // account for reservoir species
        throw std::runtime_error("Subvolume_SPDM::mapRowIndex() - invalid arguments.");
#endif
      return m_IndexerRows[i];
    };

    /**
     * Map unsorted column index to the sorted index
     * 
     * @param index Unsorted column index
     * @return Mapped column index
     */
  inline std::size_t mapColIndex(std::size_t i, std::size_t j)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i > unSpecies)||(j >= m_IndexerCols.get_cols(m_IndexerRows[i]))) // account for reservoir species
        throw std::runtime_error("Subvolume_SPDM::mapColIndex() - invalid arguments.");
#endif
      return m_IndexerCols(m_IndexerRows[i],j);
    };

    /**
     * Promote a row index in the mapping.
     * Modifies the input argument to match the resulting mapping.
     * 
     * @param index Row index
     */
  inline void moveRowUp(std::size_t & i)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i > unSpecies)) // account for reservoir species
        throw std::runtime_error("Subvolume_SPDM::moveRowUp() - invalid arguments.");
#endif
      std::swap(m_IndexerRows[i], m_IndexerRows[i - 1]);
      --i;
    };

    /**
     * Promote a column index in the mapping.
     * Modifies the input argument to match the resulting mapping.
     * 
     * @param index column index
     */
  inline void moveColLeft(std::size_t i, std::size_t & j)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i > unSpecies)||(j >= m_IndexerCols.get_cols(m_IndexerRows[i]))) // account for reservoir species
        throw std::runtime_error("Subvolume_SPDM::moveColLeft() - invalid arguments.");
#endif
      m_IndexerCols.swap(m_IndexerRows[i], j, j - 1);
      --j;
    };
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_SPDM_HPP_ */
