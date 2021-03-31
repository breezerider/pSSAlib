/**
 * @file VolumeDecomposition.hpp
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
 * Declares type for decomposition into subvolumes
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_VOLUMEDECOMPOSITION_HPP_
#define PSSALIB_DATAMODEL_DETAIL_VOLUMEDECOMPOSITION_HPP_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{

  // Forward declaration
  class DataModel;
  
namespace detail
{
  // Boundary conditions
  typedef enum tagBoundaryConditionsType
  {
    BC_Invalid,         //!<Invalid Boundary condition
    BC_Reflexive,       //!<Reflexive Boundary condition
    BC_Periodic         //! Periodic Boundary condition
  } BoundaryConditionsType;

  // Population initializer
  typedef enum tagInitialPopulationType
  {
    IP_Invalid,         //!<Invalid initial population
    IP_Distribute,      //!<Uniformly distributed initial population
    IP_Concentrate,     //!<Population concentrated in the middle
    IP_Multiply,        //!<Get the initial amount from the SBML file and uses this value for all cells
    IP_UserDefined      //!<User Defined (need callback)
  } InitialPopulationType;

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_VOLUMEDECOMPOSITION_HPP_ */

