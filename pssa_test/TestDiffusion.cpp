/**
 * @file TestDiffusion.cpp
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
 * Implementation ofr the diffusion test
 */

#include "TestDiffusion.h"

TestDiffusion::TestDiffusion()
{
  // Do nothing
}

TestDiffusion::~TestDiffusion()
{
  // Do nothing
}

bool TestDiffusion::Setup()
{
  std::string inputFile = "sbml/Diffusion.sbml";
  if (!m_SimInfo.readSBMLFile(inputFile))
  {
    std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
    return false;
  }

  m_SimInfo.eBoundaryConditions = pssalib::datamodel::detail::BC_Periodic;
  m_SimInfo.setDims(2, 3, 3);
  m_SimInfo.eInitialPopulation = pssalib::datamodel::detail::IP_Concentrate;

  return TestBase::Setup();
}

void TestDiffusion::ReactionCallback(pssalib::datamodel::DataModel * dm, REAL t)
{
  INTEGER totalA = 0, totalB = 0;
  INTEGER maxA = 0, maxB = 0;
  for (UINTEGER svi = 0; svi < dm->getSubvolumesCount(); ++svi)
  {
    pssalib::datamodel::detail::Subvolume & subVol = dm->getSubvolume(svi);
    totalA += subVol.population(1);
    totalB += subVol.population(2);
    if (subVol.population(1) > maxA) maxA = subVol.population(1);
    if (subVol.population(2) > maxB) maxB = subVol.population(2);
  }

  // Totals shouldn't change.
  if (totalA != 100 || totalB != 100)
  {
    std::cerr << "TestDiffusion::ReactionCallback: Incorrect number of molecules, A=" << totalA << " B=" << totalB << std::endl;
  }

  // By t=1, molecules should have distributed among the subvolumes.
  if (t > 1.0)
  {
    if (maxA > 50 || maxB > 50)
    {
      std::cerr << "TestDiffusion::ReactionCallback: Molecules do not seem to be diffusing, maxA=" << maxA << " maxB=" << maxB << std::endl;
    }
  }

  TestBase::ReactionCallback(dm, t);
}
