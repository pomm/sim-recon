// $Id$
//
//    File: DEventRFBunch_factory_Combo.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifndef _DEventRFBunch_factory_Combo_
#define _DEventRFBunch_factory_Combo_

#include <iostream>
#include <deque>

#include "JANA/JFactory.h"
#include "particleType.h"

#include "TRACKING/DTrackTimeBased.h"
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>

#include "PID/DEventRFBunch.h"
#include "PID/DParticleID.h"
#include "PID/DDetectorMatches.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DDetectorMatches.h"

#include "ANALYSIS/DParticleComboBlueprint.h"

using namespace jana;
using namespace std;

class DEventRFBunch_factory_Combo:public jana::JFactory<DEventRFBunch>
{
	public:
		DEventRFBunch_factory_Combo(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventRFBunch_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DParticleID* dParticleID;

		double dRFBunchFrequency;
		double dTargetCenterZ;
		double dTargetRadius;
		double dTargetLength;

		void Get_StartTime(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, double& locStartTime, double& locStartTimeVariance);
		void Calc_StartTime(const DNeutralShower* locNeutralShower, DVector3 locVertex, double& locStartTime, double& locStartTimeVariance);
		double Calc_StartTimeVariance(const DNeutralShower* locNeutralShower, const DVector3& locPathVector);
};

#endif // _DEventRFBunch_factory_Combo_

