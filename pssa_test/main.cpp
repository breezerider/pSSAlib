/**
 * @file main.cpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Aug 2015
 * @section LICENSE
 * 
 * The GPLv2 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Main test executable
 */

#include <sbml/SBMLTypes.h>
LIBSBML_CPP_NAMESPACE_USE

#include "PSSA.h"

#include "TestDiffusion.h"
#include "TestReaction.h"
#include "TestReactionDiffusion.h"

int main(int argc, char** argv)
{
  PSSALIB_MPI_IO_INIT;
  
  PSSALIB_MPI_COUT_OR_NULL << "Running pure diffusion test..." << std::endl;

  TestBase* test_diffusion = new TestDiffusion;
  if (!test_diffusion->Test())
    PSSALIB_MPI_CERR_OR_NULL << "Diffusion test failed!" << std::endl;

  PSSALIB_MPI_COUT_OR_NULL << "Running pure reaction test..." << std::endl;

  TestBase* test_reaction = new TestReaction;
  if (!test_reaction->Test())
    PSSALIB_MPI_CERR_OR_NULL << "Reaction test failed!" << std::endl;

  PSSALIB_MPI_COUT_OR_NULL << "Running reaction-diffusion test..." << std::endl;

  TestBase* test_reactiondiffusion = new TestReactionDiffusion;
  if (!test_reactiondiffusion->Test())
    PSSALIB_MPI_CERR_OR_NULL << "Reaction-diffusion test failed!" << std::endl;

  return 0;
}
