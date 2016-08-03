/*
 * DIRC_init.cc
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 *  Modified on: 0ct 7, 2013, yqiang, added DIRCTruthHit factory
 *  Modified on: Aug 2, 2016, mpatsyuk, added DIRCMirHit factory
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DDIRCHit.h"
#include "DDIRCTruthHit.h"
#include "DDIRCTruthPoint.h"
#include "DDIRCMirHit.h"

jerror_t DIRC_init(JEventLoop *loop) {
	/// Create and register DIRC data factories
	loop->AddFactory(new JFactory<DDIRCHit>());
	loop->AddFactory(new JFactory<DDIRCTruthHit>());
	loop->AddFactory(new JFactory<DDIRCTruthPoint>());
	loop->AddFactory(new JFactory<DDIRCMirHit>());

	return NOERROR;
}

