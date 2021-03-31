/**
 * @file SpeciesReference.cpp
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
 * Implements a container for certain properties of SBML SpeciesReference objects
 * that are relevant for simulations with Stochastic Simulation Algorithms
 */

#include "../../../include/datamodel/detail/SpeciesReference.h"
#include "../../../include/datamodel/detail/Reaction.h"
#include "../../../include/datamodel/detail/Species.h"
#include "../../../include/datamodel/detail/Model.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
#ifdef HAVE_LIBSBML
    /*
     * Assign properties from an SBML SpeciesReference object.
     * 
     * Implementation.
     */
    bool SpeciesReference::assign(const LIBSBML_CPP_NAMESPACE::SpeciesReference * speciesReference, Reaction * reaction)
    {
      // store the pointer to container
      m_ptrReaction = reaction;

      // call the base class method
      reinterpret_cast<Base*>(this)->assign(speciesReference);

      // validate the pointer to SBML SpeciesReference object
      if((NULL != speciesReference)&&(NULL != m_ptrReaction))
      {
        // species index
        setIndex(m_ptrReaction->getModel()->getSpeciesIndexById(speciesReference->getSpecies()));

        // species attributes
        setConstant(m_ptrReaction->getModel()->getSpecies(getIndex())->isConstant());

        // species stoichiometry
        if(round(fabs(speciesReference->getStoichiometry())) < std::numeric_limits<BYTE>::max())
          setStoichiometry((BYTE)round(fabs(speciesReference->getStoichiometry())));
        else
          return false;
      }
      else // default to reservoir species
      {
        // species index
        setIndex(0);
        // species attributes
        setConstant(false);
        // species stoichiometry
        setStoichiometry((BYTE)0);
      }

      return true;
    }
#endif
    /*
     * Get a string represantation of this object.
     * 
     * Implementation.
     */
    STRING SpeciesReference::toString() const
    {
      STRINGSTREAM ssTemp;
      ssTemp << getStoichiometryAbs() << " * ";
      if(NULL != m_ptrReaction)
        if(PSSALIB_SPECIES_ID_RESERVOIR == m_unSpeciesIndex)
          ssTemp << "[ ]";
        else
          ssTemp << m_ptrReaction->getModel()->getSpecies(getIndex())->toString();
      else
        ssTemp << Base::toString();
      return ssTemp.str();
    }

} } } // close namespaces detail, datamodel & pssalib
