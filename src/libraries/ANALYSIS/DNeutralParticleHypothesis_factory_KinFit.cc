// $Id$
//
//    File: DNeutralParticleHypothesis_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DNeutralParticleHypothesis_factory_KinFit.h"

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	locEventLoop->GetSingle(dParticleID);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DNeutralParticleHypothesis_factory_KinFit::evnt()");
#endif

 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	map<const DKinFitParticle*, DNeutralParticleHypothesis*> locKinFitParticleMap;

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		set<const DParticleCombo*> locParticleCombos;
		locKinFitResultsVector[loc_i]->Get_ParticleCombos(locParticleCombos);

		set<const DParticleCombo*>::iterator locComboIterator = locParticleCombos.begin();
		const DParticleCombo* locParticleCombo = *locComboIterator;

		map<const DKinematicData*, const DKinFitParticle*> locReverseParticleMapping;
		locKinFitResultsVector[loc_i]->Get_ReverseParticleMapping(locReverseParticleMapping);
		for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
		{
			const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
			{
				if(!locParticleComboStep->Is_FinalParticleDetected(loc_k))
					continue;
				if(!locParticleComboStep->Is_FinalParticleNeutral(loc_k))
					continue;
				const DKinFitParticle* locKinFitParticle = locReverseParticleMapping[locParticleComboStep->Get_FinalParticle_Measured(loc_k)];

				map<const DKinFitParticle*, DNeutralParticleHypothesis*>::iterator locNewHypoIterator = locKinFitParticleMap.find(locKinFitParticle);
				if(locNewHypoIterator != locKinFitParticleMap.end())
				{
					for(locComboIterator = locParticleCombos.begin(); locComboIterator != locParticleCombos.end(); ++locComboIterator)
						locNewHypoIterator->second->AddAssociatedObject(*locComboIterator);
					continue; //new particle already created for this kinfit particle
				}

				const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_k));
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticleComboStep->Get_FinalParticle(loc_k));

				DNeutralParticleHypothesis* locNewNeutralParticleHypothesis = Build_NeutralParticleHypothesis(locNeutralParticleHypothesis, locKinFitParticle, locNeutralShower, locParticleCombo);
				locKinFitParticleMap[locKinFitParticle] = locNewNeutralParticleHypothesis;

				for(locComboIterator = locParticleCombos.begin(); locComboIterator != locParticleCombos.end(); ++locComboIterator)
					locNewNeutralParticleHypothesis->AddAssociatedObject(*locComboIterator);

				_data.push_back(locNewNeutralParticleHypothesis);
			}
		}
	}

	return NOERROR;
}

DNeutralParticleHypothesis* DNeutralParticleHypothesis_factory_KinFit::Build_NeutralParticleHypothesis(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DKinFitParticle* locKinFitParticle, const DNeutralShower* locNeutralShower, const DParticleCombo* locParticleCombo)
{
	DNeutralParticleHypothesis* locNewNeutralParticleHypothesis = new DNeutralParticleHypothesis(*locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralShower);

 	vector<const JObject*> locObjects;
	locNeutralParticleHypothesis->GetT(locObjects);
	for(size_t loc_i = 0; loc_i < locObjects.size(); ++loc_i)
		locNewNeutralParticleHypothesis->AddAssociatedObject(locObjects[loc_i]);

	locNewNeutralParticleHypothesis->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewNeutralParticleHypothesis->setPosition(DVector3(locKinFitParticle->Get_CommonVertex().X(),locKinFitParticle->Get_CommonVertex().Y(),locKinFitParticle->Get_CommonVertex().Z()));
	locNewNeutralParticleHypothesis->setTime(locKinFitParticle->Get_CommonTime());
	locNewNeutralParticleHypothesis->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());

	if(locKinFitParticle->Get_ShowerEnergy() > 0.0) //particle was used in the fit as a neutral shower
		locNewNeutralParticleHypothesis->setPathLength(locKinFitParticle->Get_PathLength(), locKinFitParticle->Get_PathLengthUncertainty());
	else
	{
		double locPathLength = locNewNeutralParticleHypothesis->pathLength() - locKinFitParticle->Get_PathLength();
		double locPathLengthUncertainty_Orig = locNewNeutralParticleHypothesis->pathLength_err();
		double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
		double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
		locNewNeutralParticleHypothesis->setPathLength(locPathLength, locPathLengthUncertainty);
	}

	// Calculate DNeutralParticleHypothesis FOM
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();
	double locPropagatedRFTime = dParticleID->Calc_PropagatedRFTime(locNewNeutralParticleHypothesis, locEventRFBunch);
	locNewNeutralParticleHypothesis->setT0(locPropagatedRFTime, sqrt(locEventRFBunch->dTimeVariance), locEventRFBunch->dTimeSource);

	// Calculate DNeutralParticleHypothesis FOM
	unsigned int locNDF = 0;
	double locChiSq = 0.0;
	double locFOM = -1.0; //undefined for non-photons
	if(locNewNeutralParticleHypothesis->PID() == Gamma)
	{
		double locTimePull = 0.0;
		locChiSq = dParticleID->Calc_TimingChiSq(locNewNeutralParticleHypothesis, locNDF, locTimePull);
		locFOM = TMath::Prob(locChiSq, locNDF);
	}
	locNewNeutralParticleHypothesis->dChiSq = locChiSq;
	locNewNeutralParticleHypothesis->dNDF = locNDF;
	locNewNeutralParticleHypothesis->dFOM = locFOM;

	return locNewNeutralParticleHypothesis;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::fini(void)
{
	return NOERROR;
}


