/**
 * @file TestBase.h
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
 * base testing class
 */

#pragma once

#include "PSSA.h"

class TestBase
{
public:
	TestBase();
	virtual ~TestBase();

	bool Test();

protected:
	virtual bool Setup();
	virtual void ReactionCallback(pssalib::datamodel::DataModel* dm, REAL t);

	pssalib::datamodel::SimulationInfo m_SimInfo;
	pssalib::PSSA *pssa;

	friend void reaction_callback_wrapper(pssalib::datamodel::DataModel* dm, REAL t, void* user);
};
