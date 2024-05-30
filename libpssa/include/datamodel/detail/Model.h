/**
 * @file Model.h
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
 * Declares a container for properties of SBML Model objects that are
 * relevant for simulations using Stochastic Simulation Algorithms
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_MODEL_HPP_
#define PSSALIB_DATAMODEL_DETAIL_MODEL_HPP_

#include "Base.hpp"
#include "../../util/InplaceMemory.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  // Forward declarations
  class Species;
  class Reaction;
  class SBMLHelper;

  typedef HASHMAP <STRING, UINTEGER> MAP_ID2IDX;

  /**
   * @class Model
   * @brief Declares a container for an SBML model properties.
   */
  class Model : public Base
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagModelFlags
    {
      mfShallowCopy = (bfAll + 1),
      mfCompartmentVolumeSet = mfShallowCopy << 1,
      mfDelaysSet = mfCompartmentVolumeSet << 1,

      mfAll = mfCompartmentVolumeSet + mfDelaysSet + mfShallowCopy
    } ModelFlags;

  ////////////////////////////////
  // Attributes
  protected:
    //! Model species
    Species  *m_arSpecies;
    UINTEGER m_unSpecies;

    //! Mapping from species ids to indexes
    MAP_ID2IDX m_mapSpeciesId2Index;

    //! Model reactions
    Reaction *m_arReactions;
    UINTEGER m_unReactions,
             m_unDiffusionReactions;

    //! Compartment
    REAL     m_dCompartmentVolume;
    BYTE     m_uVolumeDims;

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    Model()
      : Base()
      , m_arSpecies(NULL)
      , m_unSpecies(0)
      , m_arReactions(NULL)
      , m_unReactions(0)
      , m_unDiffusionReactions(0)
      , m_dCompartmentVolume(0.0)
      , m_uVolumeDims(0)
    {
      // Do nothing
    };

    //! Destructor
  virtual ~Model()
    {
      // Free resources
      free();
    };

  /////////////////////////////////////
  // Methods
  public:

    /**
     * Release any allocated datastructures
     */
  virtual void free();

    /**
     * Reset all properties' values.
     */
  virtual void unset()
    {
      // call the base class method
      Base::unset();

      m_unSpecies = 0;
      m_unReactions = 0;
      m_unDiffusionReactions = 0;
      m_dCompartmentVolume = 0.0;
      m_uVolumeDims = 0;
    };
#ifdef HAVE_LIBSBML
    /**
     * Assign properties from an SBML Model object.
     * 
     * @param model an SBML @link Model object.
     */
    bool assign(const LIBSBML_CPP_NAMESPACE::Model * model, SBMLHelper & helper, STRING strCompartmentId = STRING());
#endif
    /**
     * Set compartment volume.
     * 
     * @param compartmentVolume new value.
     */
    void setCompartmentVolume(REAL compartmentVolume)
    {
      m_unFlags |= mfCompartmentVolumeSet;
      m_dCompartmentVolume = compartmentVolume;
    }

    /**
     * Get compartment volume.
     */
    REAL getCompartmentVolume() const
    {
      if(m_unFlags & mfCompartmentVolumeSet)
        return m_dCompartmentVolume;
      else
        return 0.0;
    }

    /**
     * Get compartment volume dimensions.
     */
    BYTE getCompartmentVolumeDimensions() const
    {
      if(m_unFlags & mfCompartmentVolumeSet)
        return m_uVolumeDims;
      else
        return 0;
    }

    /**
     * Set whether this model contains delayed reactions.
     * 
     * @return @true if this is a boundary condition, @false otherwise.
     */
  inline void setDelays(bool delays) { delays ? (m_unFlags |= mfDelaysSet) : (m_unFlags &= ~mfDelaysSet); };

    /**
     * Query whether this model contains delayed reactions.
     * 
     * @return @true if this model has delayed reactions, @false otherwise.
     */
  inline bool isDelaysSet() const { return (m_unFlags & mfDelaysSet); };
#ifdef HAVE_LIBSBML
    /**
     * Add a species to the model.
     * 
     * @param i A zero-based index of the species in the SBML model.
     * @param species An SBML Species object.
     * @param diffusionConstant The diffusion constant associated with this species.
     * @return @true on success, @false otherwise.
     */
    bool assignSpecies(UINTEGER i, const LIBSBML_CPP_NAMESPACE::Species * species, SBMLHelper & helper);
#endif
    /**
     * Allocate species vector for the model.
     */
    void allocSpecies(UINTEGER unSpecies);

    /**
     * Get a species from the model.
     * 
     * @param i A zero-based index of the species in this model.
     * @param species An SBML Species object.
     * @param diffusionConstant The diffusion constant associated with this species.
     * @return A pointer to species object or NULL on failure.
     */
    Species * getSpecies(UINTEGER i) const;

    /**
     * Get a species from the model.
     * 
     * @param id An identifier of the SBML Species object.
     * @return index of the Species object in the model.
     */
    UINTEGER getSpeciesIndexById(STRING id) const
    {
      MAP_ID2IDX::const_iterator it = m_mapSpeciesId2Index.find(id);
      if(m_mapSpeciesId2Index.end() == it)
      {
        STRINGSTREAM ssTemp;
        ssTemp << BOOST_CURRENT_FUNCTION << " : invalid species id '" << id << "'";
        throw std::runtime_error(ssTemp.str());
      }
      return it->second;
    }

    /**
     * Get a species from the model.
     * 
     * @param species Pointer to an instance of type Species
     *                that belongs to this model.
     * @return index of the Species object in the model.
     */
    UINTEGER getSpeciesIndex(const Species * species) const;

    /**
     * Get number of species.
     * 
     * @return Number of species.
     */
    UINTEGER getSpeciesCount() const
    {
      return m_unSpecies;
    }
#ifdef HAVE_LIBSBML
    /**
     * Add a reaction to the model.
     * 
     * @param i A zero-based index of the reaction in the SBML model.
     * @param reaction An SBML Reaction object.
     * @param diffusionConstant The diffusion constant associated with this species.
     * @return @true on success, @false otherwise.
     */
    bool assignReaction(UINTEGER i, const LIBSBML_CPP_NAMESPACE::Reaction * reaction, SBMLHelper & helper);
#endif
    /**
     * Allocate reactions vector for the model.
     */
    void allocReactions(UINTEGER unReactions);

    /**
     * Get a species from the model.
     * 
     * @param i A zero-based index of the species in this model.
     * @param species An SBML Species object.
     * @param diffusionConstant The diffusion constant associated with this species.
     * @return A pointer to species object or NULL on failure.
     */
    Reaction * getReaction(UINTEGER i) const;

    /**
     * Get number of species.
     * 
     * @return Number of species.
     */
    UINTEGER getReactionsCount() const
    {
      return m_unReactions;
    }

    /**
     * Shallow copy attribute values from another instance
     */
  virtual void copy(const Model & other)
    {
      Base::copy(other);

      m_arSpecies = other.m_arSpecies;
      m_unSpecies = other.m_unSpecies;
      m_mapSpeciesId2Index = other.m_mapSpeciesId2Index;
      m_arReactions = other.m_arReactions;
      m_unReactions = other.m_unReactions;
      m_dCompartmentVolume = other.m_dCompartmentVolume;

      m_unFlags |= mfShallowCopy;
    }

    /**
     * Swap attribute values with another instance
     */
  virtual void swap(Model & other)
    {
      Base::swap(other);

      std::swap(m_arSpecies, other.m_arSpecies);
      std::swap(m_unSpecies, other.m_unSpecies);
      std::swap(m_mapSpeciesId2Index, other.m_mapSpeciesId2Index);
      std::swap(m_arReactions, other.m_arReactions);
      std::swap(m_unReactions, other.m_unReactions);
      std::swap(m_dCompartmentVolume, other.m_dCompartmentVolume);
    }

    /**
     * Normalize reactions in the model.
     */
    void normalize();
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_MODEL_HPP_ */
