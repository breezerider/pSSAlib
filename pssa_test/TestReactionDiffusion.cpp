/**
 * @file TestReactionDiffusion.cpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Aug 2015
 * @section LICENSE
 * 
 * The GPLv2 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Implementation of the reaction-diffusion test
 */

#include "TestReactionDiffusion.h"

TestReactionDiffusion::TestReactionDiffusion()
{
}

TestReactionDiffusion::~TestReactionDiffusion()
{
}

bool TestReactionDiffusion::Setup()
{
  std::string inputFile = "sbml/Multimerization.sbml";
  if (!m_SimInfo.readSBMLFile(inputFile))
  {
    std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
    return false;
  }

  m_SimInfo.setDims(2, 3, 3);
  m_SimInfo.eBoundaryConditions = pssalib::datamodel::detail::BC_Periodic;
  m_SimInfo.eInitialPopulation = pssalib::datamodel::detail::IP_Concentrate;

  return TestBase::Setup();
}

void TestReactionDiffusion::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
  INTEGER total = 0;
  INTEGER max = 0;
  for (UINTEGER svi = 0; svi < dm->getSubvolumesCount(); ++svi)
  {
    pssalib::datamodel::detail::Subvolume & subVol = dm->getSubvolume(svi);
    INTEGER thiscount = 0;
    thiscount += subVol.population(1);
    thiscount += subVol.population(2) * 2;
    thiscount += subVol.population(3) * 3;
    thiscount += subVol.population(4) * 4;
    thiscount += subVol.population(5) * 5;

    total += thiscount;
    if (thiscount > max) max = thiscount;
  }

  // Totals shouldn't change.
  if (total != 100)
  {
    std::cerr << "TestReactionDiffusion::ReactionCallback: Incorrect number of monomers, total=" << total << std::endl;
  }

  // By t=1, molecules should have distributed among the subvolumes.
  if (t > 1.0)
  {
    if (max > 70 || max > 70)
    {
      std::cerr << "TestReactionDiffusion::ReactionCallback: Molecules do not seem to be diffusing, max=" << max << std::endl;
    }
  }

  TestBase::ReactionCallback(dm, t);
}
