/**
 * @file DataModels.cpp
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
 * Implements all the data structures required library's methods
 */

#include "../../include/PSSA.h"
#include "../../include/util/Indexing.h"
#include "../../include/util/InplaceMemory.h"

#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/DataModel_DM.h"
#include "../../include/datamodel/DataModel_PDM.h"
#include "../../include/datamodel/DataModel_SPDM.h"
#include "../../include/datamodel/DataModel_PSSACR.h"

namespace pssalib
{
  namespace datamodel
  {
    ////////////////////////////////////////
    // Generic data model class

    //! Default constructor
    DataModel::DataModel() :
        m_unReactionWrappers(0)
      , m_arReactionWrappers(NULL)
      , m_unSubvolumes(0)
      , m_arSubvolumes(NULL)
      , m_uDims(0)
      , m_arunDims(0)
      , dTotalPropensity(0.0)
      , mu(0)
      , nu(0)
      , nu_D(0)
    {
      // Do nothing
    }

    //! Destructor
    DataModel::~DataModel()
    {
      // Free resources
      free();
    }

    //! Free allocated datastructures
    void DataModel::free()
    {
      if(NULL != m_arReactionWrappers)
        util::inplace_free(m_unReactionWrappers, &m_arReactionWrappers);

      if(NULL != m_arSubvolumes)
      {
        freeSubvolumes();
        m_arSubvolumes = NULL;
        m_unSubvolumes = 0;
      }

      if(NULL != m_arunDims)
      {
        delete [] m_arunDims;
        m_arunDims = NULL;
        m_uDims = 0;
      }

      // call base class method
      detail::Model::free();
    }

    void DataModel::setup(BYTE dims, const UINTEGER *pDims, const detail::BoundaryConditionsType & bc)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      // check arguments
      if((0 != dims)&&(NULL == pDims))
        throw std::runtime_error("DataModel::setup() - invalid arguments.");
      if((0 == m_unReactions)&&(0 == m_unSpecies))
        throw std::runtime_error("DataModel::setup() - invalid model definition.");
#endif
      // pnormalize the system of reaction
      for(UINTEGER ri = 0; ri < getReactionsCount(); ++ri)
        getReaction(ri)->normalize();

      // count the number of subvolumes
      UINTEGER subvolumes = 1;
      for(BYTE di = 0; di < dims; ++di)
        subvolumes *= pDims[di];

      // check arguments
      if((dims > 0)&&(0 == subvolumes))
        throw std::runtime_error("DataModel::setup() - singleton dimensions are not allowed.");

      m_uDims = dims;
      if(pDims != m_arunDims)
      {
        if(NULL != m_arunDims)
          delete [] m_arunDims;

        if(NULL != pDims)
        {
          m_arunDims = new UINTEGER[m_uDims];
          memcpy(m_arunDims, pDims, sizeof(UINTEGER)*m_uDims);
        }
        else
          m_arunDims = NULL;
      }

      // initialize internal model representation
      setupReactionWrappers(subvolumes);

      // initialize subreaction decomposition
      setupVolumeDecomposition(subvolumes, bc);
    }

    /*
     * Initializes subvolume decomposition, i.e., neighbourhood lists, aboundary conditions
     */
    void DataModel::setupVolumeDecomposition(UINTEGER subvolumes, const detail::BoundaryConditionsType & bc)
    {
      boost::scoped_array<UINTEGER> arunSub(new UINTEGER[m_uDims]);

      // free previously allocated memory
      if(NULL != m_arSubvolumes)
        freeSubvolumes();
      m_unSubvolumes = subvolumes;
      m_arSubvolumes = new detail::Subvolume *[m_unSubvolumes];
      memset(m_arSubvolumes, 0, sizeof(detail::Subvolume *)*m_unSubvolumes);

      if(m_uDims > 0)
      {
        // Boundary conditions
        boost::scoped_ptr<tagBCHelper> bcHelper;
        switch(bc)
        {
        case detail::BC_Periodic:
          bcHelper.reset(new struct tagPeriodicBCHelper);
        break;
        case detail::BC_Reflexive:
        {
          m_unFlags &= dmfBCReflexive;
          bcHelper.reset(new struct tagReflexiveBCHelper);
        }
        break;
        default:
        {
          throw std::runtime_error("DataModel::setupVolumeDecomposition() - invalid boundary conditions.");
        }
        break;
        }

        for(UINTEGER svi = 0; svi < m_unSubvolumes; ++svi)
        {
          detail::Subvolume * pSV = allocateSubvolume();

          pSV->allocate(m_unReactionWrappers, m_unSpecies, m_uDims);

          util::ind2sub(m_uDims, m_arunDims, svi, arunSub.get());

          for(BYTE di = 0; di < m_uDims; ++di)
          {
            UINTEGER unCurrSub = arunSub[di];

            arunSub[di] = bcHelper.get()->prev(unCurrSub, m_arunDims[di]);
            util::sub2ind(m_uDims, m_arunDims, arunSub.get(), pSV->arNeighbouringSubvolumes[2*di]);

            arunSub[di] = bcHelper.get()->next(unCurrSub, m_arunDims[di]);
            util::sub2ind(m_uDims, m_arunDims, arunSub.get(), pSV->arNeighbouringSubvolumes[2*di + 1]);

            arunSub[di] = unCurrSub;
          }

          m_arSubvolumes[svi] = pSV;
        }
      }
      else
      {
        detail::Subvolume * pSV = allocateSubvolume();

        pSV->allocate(m_unReactionWrappers, m_unSpecies, 0);
        m_arSubvolumes[0] = pSV;
      }
    }

    /*
     * Sets up reaction wrappers
     */
    void DataModel::setupReactionWrappers(UINTEGER subvolumes)
    {
      // clean up
      if(NULL != m_arReactionWrappers)
        util::inplace_free(m_unReactionWrappers, &m_arReactionWrappers);

      m_unReactionWrappers = 0;
      REAL subreactorVolume = ((subvolumes > 1) ? m_dCompartmentVolume / ((REAL)subvolumes) : m_dCompartmentVolume), dH2inv = ((m_uDims > 0) ? pow(subreactorVolume, - 2.0 / (std::max((REAL)m_uDims, 2.0))) : m_dCompartmentVolume);

      // Reset minimal propensity
      crsdVolume.minValue = std::numeric_limits<REAL>::max();

      // Chemical reactions
      for(UINTEGER ri = 0, pass = 0; ri < getReactionsCount(); ++ri)
      {
        detail::Reaction * r = getReaction(ri);

        do
        {
          INTEGER exponent = 0;
          UINTEGER factor  = 1;
          for(UINTEGER rri = 0; rri < ((0 == pass) ? r->getReactantsCount() : r->getProductsCount()); ++rri)
          {
            const detail::SpeciesReference * sr = ((0 == pass) ? r->getReactantsListAt(rri) : r->getProductsListAt(rri));
            if((NULL != sr)&&(!sr->isReservoir()))
            {
              exponent -= sr->getStoichiometry();
              factor   *= gsl_sf_fact(sr->getStoichiometryAbs());
            }
          }

          REAL temp = pow(subreactorVolume, exponent+1) * ((REAL)factor);
          if(0 == pass)
          {
            temp *= r->getForwardRate();
            r->setForwardRate(temp);
          }
          else
          {
            temp *= r->getReverseRate();
            r->setReverseRate(temp);
          }

          temp /= ((REAL)factor);
          if((temp > 0.0)&&(temp < crsdVolume.minValue))
            crsdVolume.minValue = temp;

          if(!(r->isReversible()))
            break;
        }
        while(pass++ < 1);

        if(0 != pass)
          m_unReactionWrappers += 2;
        else
          ++m_unReactionWrappers;
        pass = 0;
      }

      // Species diffusion
      if(1 < subvolumes)
      {
        for(UINTEGER si = 0; si < getSpeciesCount(); ++si)
        {
          detail::Species * s = getSpecies(si);
          if(s->isSetDiffusionConstant())
          {
            REAL temp = s->getDiffusionConstant() * dH2inv;
            if (temp < crsdVolume.minValue)
              crsdVolume.minValue = temp;
            s->setDiffusionConstant(temp);
            ++m_unReactionWrappers;
          }
        }
      }

      if(0 == m_unReactionWrappers)
      {
        throw std::runtime_error("DataModel::setupReactionWrappers() - no reactions defined.");
      }

      // Allocate reaction wrappers
      m_arReactionWrappers = util::inplace_alloc<detail::ReactionWrapper>(m_unReactionWrappers);

      UINTEGER rwi = 0;
      // Chemical reactions
      for(UINTEGER ri = 0, pass = 0; ri < getReactionsCount(); ++ri)
      {
        do
        {
          new (m_arReactionWrappers + rwi) detail::ReactionWrapper(getReaction(ri), rwi, 1 == pass); ++rwi;
          if(!(getReaction(ri)->isReversible()))
            break;
        }
        while(pass++ < 1);
        pass = 0;
      }

      // Species diffusion
      if(1 < subvolumes)
      {
        for(UINTEGER si = 0; si < getSpeciesCount(); ++si)
        {
          if(getSpecies(si)->isSetDiffusionConstant())
          {
            new (m_arReactionWrappers + rwi) detail::ReactionWrapper(getSpecies(si), rwi); ++rwi;
          }
        }
      }
    }

    /*
     * Prints reaction network
     */
    std::ostream & DataModel::printReactionNetwork(std::ostream & os) const
    {
      const UINTEGER cutoff = 10, oversized = 5 * cutoff;
      bool skipped = false;

      os << "Reaction network of '" << toString() << "'" << ((m_unReactionWrappers >= oversized) ? " (reduced)" : "") << ":\n";

      os << "Volume " << m_dCompartmentVolume << "; # subreactors " << m_unSubvolumes << "\n\n";

      for(UINTEGER rwi = 0; rwi < m_unReactionWrappers; rwi++)
      {
        if((!skipped)&&(m_unReactionWrappers >= oversized)&&(rwi >= cutoff))
        {
          os << "\n. . .\n. . .\n. . .\n\n";
          rwi = m_unReactionWrappers - cutoff;
          skipped = true;
        }

        if(m_arReactionWrappers[rwi].isDiffusive())
          os << "Diffusion ";

        m_arReactionWrappers[rwi].getSymbolicRepresentation(os);
        os << '\n';
      }

      return os;
    }

    ////////////////////////////////////////
    // DM data model class

    //! Default constructor
    DataModel_DM::DataModel_DM()
    {
      // Do nothing
    }

    //! Destructor
    DataModel_DM::~DataModel_DM()
    {
      free();
    }

    ////////////////////////////////////////
    // PDM data model class

    //! Default constructor
    DataModel_PDM::DataModel_PDM()
    {
      // Do nothing
    }

    //! Destructor
    DataModel_PDM::~DataModel_PDM()
    {
      free();
    }

    ////////////////////////////////////////
    // SPDM data model class

    //! Default constructor
    DataModel_SPDM::DataModel_SPDM()
      : rowIndex(0)
      , colIndex(0)
    {
      // Do nothing
    }

    //! Destructor
    DataModel_SPDM::~DataModel_SPDM()
    {
      free();
    }

    ////////////////////////////////////////
    // PSRD-CR data model class

    //! Default constructor
    DataModel_PSSACR::DataModel_PSSACR()
    {
      // Do nothing
    }

    //! Destructor
    DataModel_PSSACR::~DataModel_PSSACR()
    {
      free();
    }
  }
}
