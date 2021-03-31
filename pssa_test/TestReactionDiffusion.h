/**
 * @file TestReactionDiffusion.h
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
 * Reaction-diffusion test
 */

#pragma once

#include "TestBase.h"

class TestReactionDiffusion : public TestBase
{
public:
  TestReactionDiffusion();
  virtual ~TestReactionDiffusion();

private:
  virtual bool Setup();
  virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);
};
