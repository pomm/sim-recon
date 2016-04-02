// $Id$
//
//    File: DBeamPhoton_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DBeamPhoton_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DBeamPhoton_factory::init(void)
{
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBeamPhoton_factory::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
    DApplication* dapp = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
    DGeometry* locGeometry = dapp->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
    dTargetCenterZ = 0.0;
    locGeometry->GetTargetZ(dTargetCenterZ);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
    DVector3 pos(0.0, 0.0, dTargetCenterZ);

    vector<const DTAGMHit*> tagm_hits;
    locEventLoop->Get(tagm_hits);

    for (unsigned int ih=0; ih < tagm_hits.size(); ++ih)
    {
        if (!tagm_hits[ih]->has_fADC) continue; // Skip TDC-only hits (i.e. hits with no ADC info.)		
        DVector3 mom(0.0, 0.0, tagm_hits[ih]->E);
        DBeamPhoton *gamma = new DBeamPhoton;
        gamma->setPID(Gamma);
        gamma->setMomentum(mom);
        gamma->setPosition(pos);
        gamma->setCharge(0);
        gamma->setMass(0);
        gamma->setTime(tagm_hits[ih]->t);
        gamma->setT0(tagm_hits[ih]->t, 0.200, SYS_TAGM);
        gamma->AddAssociatedObject(tagm_hits[ih]);
        _data.push_back(gamma);
    }

    vector<const DTAGHHit*> tagh_hits;
    locEventLoop->Get(tagh_hits);

    for (unsigned int ih=0; ih < tagh_hits.size(); ++ih)
    {
        if (!tagh_hits[ih]->has_fADC) continue; // Skip TDC-only hits (i.e. hits with no ADC info.)
        DVector3 mom(0.0, 0.0, tagh_hits[ih]->E);
        DBeamPhoton *gamma = new DBeamPhoton;
        gamma->setPID(Gamma);
        gamma->setMomentum(mom);
        gamma->setPosition(pos);
        gamma->setCharge(0);
        gamma->setMass(0);
        gamma->setTime(tagh_hits[ih]->t);
        gamma->setT0(tagh_hits[ih]->t, 0.350, SYS_TAGH);
        gamma->AddAssociatedObject(tagh_hits[ih]);
        _data.push_back(gamma);
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBeamPhoton_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamPhoton_factory::fini(void)
{
    return NOERROR;
}

