/**
 * @file TestReaction.cpp
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
 * Implementation of the reaction test
 */

#include "TestReaction.h"

TestReaction::TestReaction()
{
}

TestReaction::~TestReaction()
{
}

bool TestReaction::Setup()
{
  std::string inputFile = "sbml/Multimerization.sbml";
  if (!m_SimInfo.readSBMLFile(inputFile))
  {
    std::cerr << "Failed to load model file '" << inputFile << "'." << std::endl;
    return false;
  }

  m_SimInfo.eBoundaryConditions = pssalib::datamodel::detail::BC_Periodic;
//   simInfo.eVolumeType = pssalib::datamodel::SimulationInfo::VT_SingleSubvolume;
  // 	simInfo.dOmega = 1.0;
//   simInfo.unNumGridPoints = 1;
  m_SimInfo.eInitialPopulation = pssalib::datamodel::detail::IP_Concentrate;

	return TestBase::Setup();
}

void TestReaction::ReactionCallback(pssalib::datamodel::DataModel * dm, REAL t)
{
  INTEGER total = 0;
  for (UINTEGER svi = 0; svi < dm->getSubvolumesCount(); ++svi)
  {
    pssalib::datamodel::detail::Subvolume & subVol = dm->getSubvolume(svi);

    total += subVol.population(1);
    total += subVol.population(2) * 2;
    total += subVol.population(3) * 3;
    total += subVol.population(4) * 4;
    total += subVol.population(5) * 5;
  }

  // Totals shouldn't change.
  if (total != 100)
  {
    std::cerr << "TestReaction::ReactionCallback: Incorrect number of monomers, total=" << total << std::endl;
  }

  TestBase::ReactionCallback(dm, t);
}
