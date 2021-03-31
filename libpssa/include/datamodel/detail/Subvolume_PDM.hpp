/**
 * @file Subvolume_PDM.hpp
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
 * Declares a container for subvolume variables used by the
 * Partial Propensity Direct Method (Ramaswamy, 2009)
 */


#ifndef PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PDM_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PDM_HPP_

#include "../../stdheaders.h"
#include "JaggedMatrix.hpp"
#include "Subvolume.hpp"

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel_PDM;
  
namespace detail
{
  /**
   * @copydoc Subvolume
   */
  class Subvolume_PDM : public Subvolume
  {
  ////////////////////////////////
  // Friends
  public:
    friend class DataModel_PDM;

  ////////////////////////////////
  // Data structures
  private:

  ////////////////////////////////
  // Attributes
  protected:
    // Reactions
    //

    //! Propensity of each group
    REAL                       *m_ardLambda;
    //! Total propensity of each group
    REAL                       *m_ardSigma;

  ////////////////////////////////
  // Attributes
  public:
    // Reactions
    //

    //! Array of arrays of reaction partial propensities.
    JaggedMatrix<REAL>         arPi;

  ////////////////////////////////
  // Constructors
  public:
    //! Constructor
    Subvolume_PDM()
      : m_ardLambda(NULL)
      , m_ardSigma(NULL)
    {
      // Do nothing
    }

    //! Destructor
  virtual ~Subvolume_PDM()
    {
      // Clean-up
      free_PDM();
    }

  ////////////////////////////////
  // Methods
  private:

    /*
     * Free memory
     */
    void free_PDM() 
    {
      if(NULL != m_ardLambda)
      {
        delete [] m_ardLambda;
        m_ardLambda = NULL;
      }
      if(NULL != m_ardSigma)
      {
        delete [] m_ardSigma;
        m_ardSigma = NULL;
      }
    };

  ////////////////////////////////
  // Methods
  protected:

    /**
     * @copydoc Subvolume::free()
     */
  virtual void free()
    {
      // free memory
      free_PDM();

      // call base class method
      Subvolume::free();
    };

    /**
     * @copydoc Subvolume::allocate(UINTEGER,UINTEGER,BYTE)
     */
  virtual void allocate(UINTEGER reactions, UINTEGER species, BYTE dims)
    {
      // call base class method
      Subvolume::allocate(reactions, species, dims);

      // allocate memory
      m_ardLambda = new REAL[species + 1];
      memset(m_ardLambda, 0, sizeof(REAL)*(species + 1));
      m_ardSigma = new REAL[species + 1];
      memset(m_ardSigma, 0, sizeof(REAL)*(species + 1));
      arPi.reserve(species + 1, std::max(reactions / species, (UINTEGER)1));
    };

    /**
     * @copydoc Subvolume::clear(UINTEGER,UINTEGER)
     */
  virtual void clear(UINTEGER reactions, UINTEGER species)
    {
      memset(m_ardLambda, (unsigned char)0, sizeof(REAL)*(species + 1));
      memset(m_ardSigma, (unsigned char)0, sizeof(REAL)*(species + 1));
      arPi.clear();

      // call base class method
      Subvolume::clear(reactions, species);
    };


  ////////////////////////////////
  // Methods
  public:

    /**
     * Get Lambda value (see Ramaswamy, 2009) for a species
     * 
     * @param index Species index in the model
     * @return Lambda value
     */
  inline REAL & lambda(UINTEGER index)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= (unSpecies + 1))
        throw std::runtime_error("Subvolume_PDM::lambda() - invalid arguments.");
#endif
      return m_ardLambda[index];
    }

    /**
     * Get Sigma value (see Ramaswamy, 2009) for a species
     * 
     * @param index Species index in the model
     * @return Sigma value
     */
  inline REAL & sigma(UINTEGER index)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(index >= (unSpecies + 1))
        throw std::runtime_error("Subvolume_PDM::sigma() - invalid arguments.");
#endif
      return m_ardSigma[index];
    }
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SUBVOLUME_PDM_HPP_ */
