/**
 * @file TestBase.cpp
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
 * Implementation for the base testing class
 */

#include "TestBase.h"

void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user)
{
  TestBase* base = reinterpret_cast<TestBase*>(user);
  base->ReactionCallback(dm, t);
}

void progress_callback(UINTEGER a, UINTEGER b, SHORT c, void* user)
{
  std::cout << "\rProgress: sample " << a << "/" << b << " " << c << "% done..." << std::flush;
}

TestBase::TestBase()
  : pssa(NULL)
{
}

TestBase::~TestBase()
{
  if (pssa) {
    delete pssa;
  }
}

bool TestBase::Test()
{
  if (!Setup()) {
    return false;
  }

  bool result = pssa->run(&m_SimInfo);
  std::cout << std::endl;
  return result;
}

bool TestBase::Setup()
{
  m_SimInfo.dTimeStart = 0.0;
  m_SimInfo.dTimeStep = 0.1;
  m_SimInfo.dTimeEnd = 1000.0;
  m_SimInfo.strOutput = STRING();
  m_SimInfo.unSamplesTotal =
#ifdef HAVE_MPI
                              pssalib::getMPIWrapperInstance().getPoolSize() + 1;
#else
                              1;
#endif
  m_SimInfo.pArSpeciesIds = new std::vector<STRING>();
  m_SimInfo.unOutputFlags = ~(
    pssalib::datamodel::SimulationInfo::ofTrace      |
    pssalib::datamodel::SimulationInfo::ofFinalPops  |
    pssalib::datamodel::SimulationInfo::ofTrajectory |
    pssalib::datamodel::SimulationInfo::ofTiming     |
    pssalib::datamodel::SimulationInfo::ofTimePoints |
    pssalib::datamodel::SimulationInfo::ofStatus);
  m_SimInfo.unOutputFlags |= pssalib::datamodel::SimulationInfo::ofLog
    | pssalib::datamodel::SimulationInfo::ofInfo
    | pssalib::datamodel::SimulationInfo::ofError;
    
  m_SimInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());

  pssa = new pssalib::PSSA();
  if (!pssa)
  {
    std::cerr << "Failed to allocate the simulation engine." << std::endl;
    return false;
  }

  if(!pssa->setMethod(pssalib::PSSA::M_SPDM))
  {
    std::cerr << "Failed to set simulation method." << std::endl;
    return false;
  }

  pssa->SetReactionCallback(&reaction_callback_wrapper, this);
  pssa->SetProgressCallback(&progress_callback, this);

  return true;
}

void TestBase::ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t)
{
}
