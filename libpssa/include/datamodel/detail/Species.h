/**
 * @file Species.h
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
 * Declares a container for properties of SBML Species objects that are
 * relevant for simulations using Stochastic Simulation Algorithms
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SPECIES_H_
#define PSSALIB_DATAMODEL_DETAIL_SPECIES_H_

#include "Base.hpp"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  // Forward declarations
  class Model;
  class SBMLHelper;

  /**
   * @class Species
   * @brief Declares a datatype for chemical species
   */
  class Species : public Base
  {
  /////////////////////////////////////
  // Data structures
  public:
    //! Flags
    enum tagSpeciesFlags
    {
      sfInitialAmountSet = (bfAll + 1),

      sfConstant = sfInitialAmountSet << 1,
      sfBoundaryCondition = sfInitialAmountSet << 2,

      sfDiffusive = sfInitialAmountSet << 3,

      sfAll = sfInitialAmountSet + sfConstant + sfBoundaryCondition + sfDiffusive
    } SpeciesFlags;

  ////////////////////////////////
  // Attributes
  protected:
    //! Initial amount of species: equals to 'initialAmount'
    //! property (if set) or is  using computed as product of
    //! 'initialConcentration' property and reactor volume
    //! if both are set. Otherwise, it is set to 0.
    UINTEGER m_unInitialAmount;

//     //! Species population
//     UINTEGER *arunPopulation;

    //! Diffusion constant value
    REAL m_dDiffusionConstant;

    //! Container
    Model const & m_refModel;

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    Species(Model const & model)
      : Base()
      , m_unInitialAmount(0)
      , m_dDiffusionConstant(std::numeric_limits<REAL>::min())
      , m_refModel(model)
    {
      // Do nothing
    };

    //! Copy constructor
    Species(const Species & right)
      : Base(right)
      , m_unInitialAmount(right.m_unInitialAmount)
      , m_dDiffusionConstant(right.m_dDiffusionConstant)
      , m_refModel(right.m_refModel)
    {
      // Do nothing
    };

    //! Destructor
  virtual ~Species()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  private:
#ifdef HAVE_LIBSBML
    /**
     * Acquire propertie values from the annotation of an SBML Reaction object.
     * 
     * @param species an SBML Species object.
     */
    bool parseAnnotation(const LIBSBML_CPP_NAMESPACE::Species * species, SBMLHelper & helper);
#endif
  public:

    /**
     * Reset all properties' values.
     */
  virtual void unset()
    {
      // call the base class method
      Base::unset();
      m_unInitialAmount = 0;
      m_dDiffusionConstant = std::numeric_limits<REAL>::min();
    };
#ifdef HAVE_LIBSBML
    /**
     * Assign properties from an SBML Species object.
     * 
     * @param species an SBML Species object.
     * @param model the parent Model object.
     * @return @true if assignment was successful, @false otherwise.
     */
    bool assign(const LIBSBML_CPP_NAMESPACE::Species * species, SBMLHelper & helper);
#endif
    /**
     * Get the index of this species in the model.
     * 
     * @return current model species index.
     */
    UINTEGER getIndex() const;

    /**
     * Set the initial amount of species.
     * 
     * @param initialAmount new value.
     */
  inline void setInitialAmount(UINTEGER initialAmount)
    {
      m_unFlags |= sfInitialAmountSet;
      m_unInitialAmount = initialAmount;
    };

    /**
     * Get the diffusion constant of species.
     * 
     * @return current value.
     */
  inline UINTEGER getInitialAmount() const
    {
      if(m_unFlags & sfInitialAmountSet)
        return m_unInitialAmount;
      return 0;
    };

    /**
     * Get the diffusion constant of species.
     * 
     * @return current value.
     */
  inline REAL getDiffusionConstant() const
    {
      if(m_unFlags & sfDiffusive)
        return m_dDiffusionConstant;
      return 0.0;
    };

    /**
     * Set the diffusion constant of species.
     * 
     * @param diffusionConstant new value.
     */
  inline void setDiffusionConstant(const REAL diffusionConstant)
    {
      if(diffusionConstant > 0.0)
      {
        m_unFlags |= sfDiffusive;
        m_dDiffusionConstant = diffusionConstant;
      }
      else
      {
        m_unFlags &= ~sfDiffusive;
        m_dDiffusionConstant = std::numeric_limits<REAL>::min();
      }
    };

    /**
     * Check if the diffusion constant of species is set.
     * 
     * @return @true if property is set, @false otherwise.
     */
  inline bool isSetDiffusionConstant() const { return (bool)(m_unFlags & sfDiffusive); };

    /**
     * Set whether this species is should be kept constant.
     * 
     * @param c @true if this species is to be kept constant, @false otherwise.
     */
  inline void setConstant(bool c) { c ? (m_unFlags |= sfConstant) : (m_unFlags &= ~sfConstant); };

    /**
     * Set whether this species is a boundary condition (reservoir species).
     * 
     * @param bc @true if this is a boundary condition, @false otherwise.
     */
  inline void setBoundaryCondition(bool bc) { bc ? (m_unFlags |= sfBoundaryCondition) : (m_unFlags &= ~sfBoundaryCondition); };

    /**
     * Query whether this species is should be kept constant.
     * 
     * @return @true if this species is to be kept constant, @false otherwise.
     */
  inline bool isConstant() const { return (m_unFlags & sfConstant); };

    /**
     * Query whether this species is a boundary condition (reservoir species).
     * 
     * @return @true if this is a boundary condition, @false otherwise.
     */
  inline bool isBoundaryCondition() const { return (m_unFlags & sfBoundaryCondition); };
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SPECIES_H_ */
