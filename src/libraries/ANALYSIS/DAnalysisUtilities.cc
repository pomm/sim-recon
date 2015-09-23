#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DAnalysisUtilities.h"

DAnalysisUtilities::DAnalysisUtilities(JEventLoop* locEventLoop)
{
  // Get the particle ID algorithms
	vector<const DParticleID*> locPIDAlgorithms;
	locEventLoop->Get(locPIDAlgorithms);
	if(locPIDAlgorithms.size() < 1)
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
	dPIDAlgorithm = locPIDAlgorithms[0];

	dTargetZCenter = 65.0;
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber()) : NULL;
	if(locGeometry != NULL)
		locGeometry->GetTargetZ(dTargetZCenter);

	//For "Unused" tracks/showers
	//BEWARE: IF THIS IS CHANGED, CHANGE IN THE BLUEPRINT FACTORY AND THE EVENT WRITER ALSO!!
	dTrackSelectionTag = "PreSelect";
	dShowerSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);
}

bool DAnalysisUtilities::Check_IsBDTSignalEvent(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag, bool locIncludeDecayingToReactionFlag) const
{
#ifdef VTRACE
	VT_TRACER("DAnalysisUtilities::Check_IsBDTSignalEvent()");
#endif

	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//if locIncludeDecayingToReactionFlag = true, will test whether the thrown reaction could decay to the DReaction
		//Note that resonances, phi's, and omega's are automatically decayed
			//e.g. if DReaction or thrown is g, p -> pi+, pi-, omega, p; will instead treat it as g, p -> 2pi+, 2pi-, pi0, p (or whatever the omega decay products are)
	//e.g. g, p -> pi+, pi0, K0, Lambda can decay to g, p -> 2pi+, 2pi-, pi0, p
		//if locIncludeDecayingToReactionFlag = true, then it would be included as "Signal," if false, then background
		//locIncludeDecayingToReactionFlag should be true UNLESS you are explicitly checking all possible reactions that could decay to your channel in your BDT
			//e.g. could kinfit to g, p -> pi+, pi0, K0, Lambda and include it as a BDT variable

	DParticleCombo_factory_Thrown* dThrownComboFactory = static_cast<DParticleCombo_factory_Thrown*>(locEventLoop->GetFactory("DParticleCombo", "Thrown"));
	DReaction_factory_Thrown* dThrownReactionFactory = static_cast<DReaction_factory_Thrown*>(locEventLoop->GetFactory("DReaction", "Thrown"));

	//Replace omega & phi in DReaction with their decay products
	size_t locStepIndex = 0;
	vector<DReactionStep> locNewReactionSteps;
	const DReaction* locCurrentReaction = locReaction;
	DReaction locNewReaction("Interim"); //if needed

	do
	{
		const DReactionStep* locReactionStep = locCurrentReaction->Get_ReactionStep(locStepIndex);
		Particle_t locInitialPID = locReactionStep->Get_InitialParticleID();
		if((locInitialPID != omega) && (locInitialPID != phiMeson))
		{
			++locStepIndex;
			continue;
		}

		//is decaying phi or omega, replace it with its decay products
		deque<Particle_t> locDecayProducts;
		locReactionStep->Get_FinalParticleIDs(locDecayProducts);

		//find the production step
		for(size_t loc_i = 0; loc_i < locStepIndex; ++loc_i)
		{
			const DReactionStep* locProductionStep = locCurrentReaction->Get_ReactionStep(loc_i);
			deque<Particle_t> locPIDs;
			locProductionStep->Get_FinalParticleIDs(locPIDs);

			//search for the decaying particle. when found, replace it with its decay products
			bool locFoundFlag = false;
			for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
			{
				Particle_t locFinalStatePID = locPIDs[loc_j];
				if(locFinalStatePID != locInitialPID)
					continue;
				locPIDs.erase(locPIDs.begin() + loc_j);
				locPIDs.insert(locPIDs.begin(), locDecayProducts.begin(), locDecayProducts.end());
				locFoundFlag = true;
				break;
			}
			if(!locFoundFlag)
				continue;

			//make a new reaction step
			DReactionStep locNewReactionStep;
			locNewReactionStep.Set_InitialParticleID(locProductionStep->Get_InitialParticleID());
			locNewReactionStep.Set_TargetParticleID(locProductionStep->Get_TargetParticleID());
			int locMissingParticleIndex = locProductionStep->Get_MissingParticleIndex();
			for(int loc_j = 0; loc_j < int(locPIDs.size()); ++loc_j)
				locNewReactionStep.Add_FinalParticleID(locPIDs[loc_j], (loc_j == locMissingParticleIndex));
			locNewReactionSteps.push_back(locNewReactionStep);

			//make a new reaction
			DReaction locReactionBuild("Build");
			locReactionBuild.Clear_ReactionSteps();
			for(size_t loc_j = 0; loc_j < locCurrentReaction->Get_NumReactionSteps(); ++loc_j)
			{
				if(loc_j == loc_i)
					locReactionBuild.Add_ReactionStep(&locNewReactionSteps.back());
				else if(loc_j == locStepIndex)
					continue;
				else
					locReactionBuild.Add_ReactionStep(locCurrentReaction->Get_ReactionStep(loc_j));
			}
			locNewReaction = locReactionBuild;
			locCurrentReaction = &locNewReaction;
			break;
		}
	}
	while(locStepIndex < locCurrentReaction->Get_NumReactionSteps());

	//Get thrown steps
	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	Get_ThrownParticleSteps(locEventLoop, locThrownSteps);

	if(locThrownSteps.size() == 1) //nothing to replace: it either works or it doesn't
		return Check_ThrownsMatchReaction(locEventLoop, locCurrentReaction, locExclusiveMatchFlag);

	//build maps of counts of the thrown & DReaction decaying particles
	map<Particle_t, size_t> locNumDecayingParticles_Thrown;
	for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
	{
		if(locThrownSteps[loc_i].first == NULL)
			continue; //beam particle
		Particle_t locPID = locThrownSteps[loc_i].first->PID();
		if(locNumDecayingParticles_Thrown.find(locPID) == locNumDecayingParticles_Thrown.end())
			locNumDecayingParticles_Thrown[locPID] = 1;
		else
			++locNumDecayingParticles_Thrown[locPID];
	}

	map<Particle_t, size_t> locNumDecayingParticles_Reaction;
	for(size_t loc_i = 0; loc_i < locCurrentReaction->Get_NumReactionSteps(); ++loc_i)
	{
		Particle_t locPID = locCurrentReaction->Get_ReactionStep(loc_i)->Get_InitialParticleID();
		if(locPID == Gamma)
			continue; //beam particle
		if(locNumDecayingParticles_Reaction.find(locPID) == locNumDecayingParticles_Reaction.end())
			locNumDecayingParticles_Reaction[locPID] = 1;
		else
			++locNumDecayingParticles_Reaction[locPID];
	}

	//if there is a particle type in the DReaction that is not present in the thrown, it's gonna fail no matter what: check now
	map<Particle_t, size_t>::iterator locIterator = locNumDecayingParticles_Reaction.begin();
	for(; locIterator != locNumDecayingParticles_Reaction.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		if(locNumDecayingParticles_Thrown.find(locPID) == locNumDecayingParticles_Thrown.end())
			return false;
	}

	//loop through thrown steps: replace phi's, omega's with their decay products
		//also replace particles not in the DReaction with their decay products IF locIncludeDecayingToReactionFlag = true
	locStepIndex = 1;
	do
	{
		pair<const DMCThrown*, deque<const DMCThrown*> > locStepPair = locThrownSteps[locStepIndex];

		//check to see if the thrown decaying particle is present in the dreaction
		const DMCThrown* locThrownParent = locStepPair.first;
		Particle_t locInitialPID = locThrownParent->PID();
		bool locInitialPIDFoundFlag = (locNumDecayingParticles_Reaction.find(locInitialPID) != locNumDecayingParticles_Reaction.end());

		//if it was not found, the only way it can be a match is if the decay products ARE in the reaction step: decay it in place
			//also: omega and phi have non-negligible width: cannot constrain their peaks in kinfit or BDT: replace with their decay products
		if(((!locInitialPIDFoundFlag) && locIncludeDecayingToReactionFlag) || (locInitialPID == omega) || (locInitialPID == phiMeson))
		{
			Replace_DecayingParticleWithProducts(locThrownSteps, locStepIndex);
			if(locNumDecayingParticles_Thrown[locInitialPID] == 1)
				locNumDecayingParticles_Thrown.erase(locInitialPID);
			else
				--locNumDecayingParticles_Thrown[locInitialPID];
		}
		else
			++locStepIndex;
	}
	while(locStepIndex < locThrownSteps.size());

	if(!locIncludeDecayingToReactionFlag)
	{
		//don't try decaying thrown particles: compare as-is
		const DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locThrownSteps);
		const DParticleCombo* locThrownCombo = dThrownComboFactory->Build_ThrownCombo(locEventLoop, locThrownReaction, locThrownSteps);
		bool locCheckResult = Check_ThrownsMatchReaction(locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

		delete locThrownCombo;
		delete locThrownReaction;
		return locCheckResult;
	}

	//ok, if there are still an unequal # of parents for a given PID between thrown & DReaction (i.e. both #'s are non-zero):
		//it's too confusing to figure out which should be replaced by their decay products and which shouldn't
		//so, try all possibilities, and see if any of them match. 

	//build PIDs-to-replace vectors //easiest to manipulate a 1D vector
	vector<Particle_t> locPIDVector;
	vector<int> locResumeAtIndex;
	for(locIterator = locNumDecayingParticles_Reaction.begin(); locIterator != locNumDecayingParticles_Reaction.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		size_t locNumThrown = locNumDecayingParticles_Thrown[locPID];
		if(locNumThrown < locIterator->second)
			return false; //thrown doesn't match
		if(locNumThrown == locIterator->second)
			continue; //all is well
		for(size_t loc_i = 0; loc_i < locNumThrown - locIterator->second; ++loc_i)
		{
			locPIDVector.push_back(locPID);
			locResumeAtIndex.push_back(0);
		}
	}

	//if no additional replacements to make: check it
	if(locPIDVector.empty())
	{
		const DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locThrownSteps);
		const DParticleCombo* locThrownCombo = dThrownComboFactory->Build_ThrownCombo(locEventLoop, locThrownReaction, locThrownSteps);
		bool locCheckResult = Check_ThrownsMatchReaction(locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

		delete locThrownCombo;
		delete locThrownReaction;
		return locCheckResult;
	}

	//loop through, making all combos of particles, (almost) just like in DParticleComboBlueprint_factory
		//unlike that factory, the below may try a given combo twice: too hard to code against it though 
			//for a given PID, hard to compare current locResumeAtIndex to previous ones, since the thrown steps are continually replaced
	deque<int> locDecayReplacementIndices;
	int locParticleIndex = 0;
	//below: contains "locThrownSteps" after each replacement
		//1st deque index is replacement index, 
	deque<deque<pair<const DMCThrown*, deque<const DMCThrown*> > > > locReplacementThrownSteps(1, locThrownSteps);
	do
	{
		if(locParticleIndex == -1)
			break; //no success

		deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locCurrentThrownSteps = locReplacementThrownSteps.back();
		if(locParticleIndex == int(locPIDVector.size()))
		{
			//combo defined: try it
			const DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locCurrentThrownSteps);
			const DParticleCombo* locThrownCombo = dThrownComboFactory->Build_ThrownCombo(locEventLoop, locThrownReaction, locCurrentThrownSteps);
			bool locCheckResult = Check_ThrownsMatchReaction(locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

			delete locThrownCombo;
			delete locThrownReaction;

			if(locCheckResult)
				return true; //it worked!
			locResumeAtIndex[locParticleIndex] = 0;
			locReplacementThrownSteps.pop_back();
			--locParticleIndex;
			continue;
		}

		Particle_t locPID = locPIDVector[locParticleIndex];
		//find the next instance of this step & replace it
		bool locFoundFlag = false;
		for(size_t loc_i = locResumeAtIndex[locParticleIndex]; loc_i < locCurrentThrownSteps.size(); ++loc_i)
		{
			if(locCurrentThrownSteps[loc_i].first == NULL)
				continue;
			Particle_t locInitialPID = locCurrentThrownSteps[loc_i].first->PID();
			if(locInitialPID != locPID)
				continue;
			locFoundFlag = true;

			Replace_DecayingParticleWithProducts(locCurrentThrownSteps, loc_i);
			locReplacementThrownSteps.push_back(locCurrentThrownSteps);
			locResumeAtIndex[locParticleIndex] = loc_i + 1;

			break;
		}

		if(locFoundFlag)
			++locParticleIndex;
		else
		{
			locResumeAtIndex[locParticleIndex] = 0;
			locReplacementThrownSteps.pop_back();
			--locParticleIndex;
		}
	}
	while(true);

	return false;
}

void DAnalysisUtilities::Replace_DecayingParticleWithProducts(deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps, size_t locStepIndex) const
{
	if(locStepIndex >= locThrownSteps.size())
		return;

	//find the step where this particle is produced at
	int locInitialParticleDecayFromStepIndex = -1;
	const DMCThrown* locThrownParent = locThrownSteps[locStepIndex].first;
	if(locThrownParent == NULL)
		return;

	for(size_t loc_j = 0; loc_j < locStepIndex; ++loc_j)
	{
		for(size_t loc_k = 0; loc_k < locThrownSteps[loc_j].second.size(); ++loc_k)
		{
			if(locThrownParent != locThrownSteps[loc_j].second[loc_k])
				continue;
			locInitialParticleDecayFromStepIndex = loc_j;
			break;
		}
		if(locInitialParticleDecayFromStepIndex != -1)
			break;
	}

	if(locInitialParticleDecayFromStepIndex == -1)
	{
		cout << "ERROR: SOMETHING IS WRONG WITH DAnalysisUtilities::Replace_DecayingParticleWithProducts(). ABORTING." << endl;
		abort();
	}

	//insert the decay products where the production occured
	deque<const DMCThrown*>& locProductionStepFinalParticles = locThrownSteps[locInitialParticleDecayFromStepIndex].second;
	deque<const DMCThrown*>& locDecayStepFinalParticles = locThrownSteps[locStepIndex].second;
	for(size_t loc_j = 0; loc_j < locProductionStepFinalParticles.size(); ++loc_j)
	{
		if(locProductionStepFinalParticles[loc_j] != locThrownParent)
			continue;
		locProductionStepFinalParticles.erase(locProductionStepFinalParticles.begin() + loc_j);
		locProductionStepFinalParticles.insert(locProductionStepFinalParticles.end(), locDecayStepFinalParticles.begin(), locDecayStepFinalParticles.end());
		break;
	}

	locThrownSteps.erase(locThrownSteps.begin() + locStepIndex);
}

bool DAnalysisUtilities::Check_ThrownsMatchReaction(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag) const
{
	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//note, if you decay a final state particle (e.g. k+, pi+) in your input DReaction*, a match will NOT be found: the thrown reaction/combo is truncated
	//if locExclusiveMatchFlag = false, then allow the input DReaction to be a subset of the thrown
	const DParticleCombo* locThrownCombo = NULL;
	locEventLoop->GetSingle(locThrownCombo, "Thrown");

	return Check_ThrownsMatchReaction(locThrownCombo, locReaction, locExclusiveMatchFlag);
}

bool DAnalysisUtilities::Check_ThrownsMatchReaction(const DParticleCombo* locThrownCombo, const DReaction* locReaction, bool locExclusiveMatchFlag) const
{
#ifdef VTRACE
	VT_TRACER("DAnalysisUtilities::Check_ThrownsMatchReaction()");
#endif

	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//note, if you decay a final state particle (e.g. k+, pi+) in your input DReaction*, a match will NOT be found: the thrown reaction/combo is truncated
	//if locExclusiveMatchFlag = false, then allow the input DReaction to be a subset of the thrown

	if(locThrownCombo == NULL)
		return false;
	const DReaction* locThrownReaction = locThrownCombo->Get_Reaction();

	if(locExclusiveMatchFlag)
	{
		if(locReaction->Get_NumReactionSteps() != locThrownCombo->Get_NumParticleComboSteps())
			return false;
	}

	//build map of InitialParticleDecayFromStepIndex's for input reaction: assume that if it matters, the user wanted them in the input order
	map<size_t, int> locReactionInitialParticleDecayFromStepIndexMap; //first is step index (of locReaction) where particle is initial, second is where is final state particle
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		if(loc_i == 0)
		{
			if(locReactionStep->Get_InitialParticleID() == Gamma)
				locReactionInitialParticleDecayFromStepIndexMap[0] = -1;
		}
		//loop over final state particles, and if decaying, find the step they are a parent in
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			Particle_t locFinalStatePID = locReactionStep->Get_FinalParticleID(loc_j);
			//see if this pid is a parent in a future step
			for(size_t loc_k = loc_i; loc_k < locReaction->Get_NumReactionSteps(); ++loc_k)
			{
				if(locReaction->Get_ReactionStep(loc_k)->Get_InitialParticleID() != locFinalStatePID)
					continue;
				if(locReactionInitialParticleDecayFromStepIndexMap.find(loc_k) != locReactionInitialParticleDecayFromStepIndexMap.end())
					continue; //this step already accounted for
				locReactionInitialParticleDecayFromStepIndexMap[loc_k] = loc_i;
				break;
			}
		}
	}

	//since throwns are more organized, loop through those and back-check to dreaction
	set<size_t> locMatchedInputStepIndices; //step indices in input locReaction that are already matched-to
	map<int, int> locStepMatching; //locReaction map to thrown step map
	map<int, int> locReverseStepMatching; //thrown map to locReaction step map
	for(size_t loc_i = 0; loc_i < locThrownCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locThrownParticleComboStep = locThrownCombo->Get_ParticleComboStep(loc_i);
		int locInitialParticleDecayFromStepIndex = locThrownParticleComboStep->Get_InitialParticleDecayFromStepIndex();

		//find where to start the search for this step in locReaction
		size_t locStartSearchIndex = 0; //locReaction could be a subset of the total; start at the beginning unless ...:
		if(locReverseStepMatching.find(locInitialParticleDecayFromStepIndex) != locReverseStepMatching.end())
			locStartSearchIndex = locReverseStepMatching[locInitialParticleDecayFromStepIndex] + 1; //parent step was matched, start search for this one after it

		//loop through locReaction and try to find this thrown step
		bool locMatchFoundFlag = false;
		int locPossibleMatchIndex = -1;
		for(size_t loc_j = locStartSearchIndex; loc_j < locReaction->Get_NumReactionSteps(); ++loc_j)
		{
			const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_j);
			if(locMatchedInputStepIndices.find(loc_j) != locMatchedInputStepIndices.end())
				continue; //this step was already accounted for
			//when not exact match: allow user step to have a missing unknown particle
			if(!locThrownReaction->Get_ReactionStep(loc_i)->Are_ParticlesIdentical(locReactionStep, !locExclusiveMatchFlag))
				continue; //particles aren't the same

			//ok, now check to make sure that the parent particle in this step was produced the same way in both thrown & locReaction
			if(locReactionInitialParticleDecayFromStepIndexMap.find(loc_j) == locReactionInitialParticleDecayFromStepIndexMap.end())
			{
				//this (loc_j) parent particle's production step wasn't listed in the locReaction: locReaction is probably a subset of the total
				//a match is possible but not certain: (e.g. locReaction is pi0, eta -> pi0, pi0 and this step (loc_i) is a pi0)
				locPossibleMatchIndex = loc_j;
				continue; //keep searching in case a different one should be used instead //will use this if another isn't found
			}

			int locReactionInitialParticleDecayFromStepIndex = locReactionInitialParticleDecayFromStepIndexMap[loc_j];
			if(locReactionInitialParticleDecayFromStepIndex != -1) //locReaction is not beam particle
			{
				if(locStepMatching.find(locReactionInitialParticleDecayFromStepIndex) == locStepMatching.end())
					continue; //the step in locReaction where this (loc_j) parent particle was produced was not mapped to the thrown steps yet: this is not the step we want
				int locReactionInitialParticleDecayFromStepIndexMappedBackToThrown = locStepMatching[locReactionInitialParticleDecayFromStepIndex];
				if(locInitialParticleDecayFromStepIndex != locReactionInitialParticleDecayFromStepIndexMappedBackToThrown)
					continue; //the decaying parent particle in this step (loc_j) comes from a different step in thrown (locInitialParticleDecayFromStepIndex)/reaction: continue
			}
			else if(locInitialParticleDecayFromStepIndex != -1)
				continue; //reaction is beam but thrown is not beam

			//finally, a match! register it
			locMatchedInputStepIndices.insert(loc_j);
			locStepMatching[loc_j] = loc_i;
			locReverseStepMatching[loc_i] = loc_j;
			locMatchFoundFlag = true;
			break;
		}
		if((!locMatchFoundFlag) && (locPossibleMatchIndex != -1))
		{
			//need to use the possible match
			locMatchedInputStepIndices.insert(locPossibleMatchIndex);
			locStepMatching[locPossibleMatchIndex] = loc_i;
			locReverseStepMatching[loc_i] = locPossibleMatchIndex;
			locMatchFoundFlag = true;
		}
		if(locExclusiveMatchFlag && (!locMatchFoundFlag))
			return false; //needed an exact match and it wasn't found: bail
	}
	if(locExclusiveMatchFlag)
		return true;

	//locReaction could be a subset of thrown: check if all locReaction steps found
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locMatchedInputStepIndices.find(loc_i) == locMatchedInputStepIndices.end())
			return false; //one of the input steps wasn't matched: abort!
	}

	return true;
}

void DAnalysisUtilities::Get_UnusedChargedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DChargedTrack*>& locUnusedChargedTracks) const
{
	locUnusedChargedTracks.clear();
	locEventLoop->Get(locUnusedChargedTracks, dTrackSelectionTag.c_str());

	deque<const DChargedTrack*> locSourceChargedTracks;
	locParticleCombo->Get_DetectedFinalChargedParticles_SourceObjects(locSourceChargedTracks);

	for(size_t loc_i = 0; loc_i < locSourceChargedTracks.size(); ++loc_i)
	{
		vector<const DChargedTrack*>::iterator locIterator;
		for(locIterator = locUnusedChargedTracks.begin(); locIterator != locUnusedChargedTracks.end();)
		{
			if(locSourceChargedTracks[loc_i] != *locIterator)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedChargedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedTimeBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackTimeBased*>& locUnusedTimeBasedTracks) const
{
	locUnusedTimeBasedTracks.clear();
	locEventLoop->Get(locUnusedTimeBasedTracks);

	deque<const DChargedTrack*> locSourceChargedTracks;
	locParticleCombo->Get_DetectedFinalChargedParticles_SourceObjects(locSourceChargedTracks);

	for(size_t loc_i = 0; loc_i < locSourceChargedTracks.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locSourceChargedTracks[loc_i]->dChargedTrackHypotheses[0];

		vector<const DTrackTimeBased*>::iterator locIterator;
		for(locIterator = locUnusedTimeBasedTracks.begin(); locIterator != locUnusedTimeBasedTracks.end();)
		{
			if(locChargedTrackHypothesis->candidateid != (*locIterator)->candidateid)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedTimeBasedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedWireBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackWireBased*>& locUnusedWireBasedTracks) const
{
	locUnusedWireBasedTracks.clear();
	locEventLoop->Get(locUnusedWireBasedTracks);

	deque<const DChargedTrack*> locSourceChargedTracks;
	locParticleCombo->Get_DetectedFinalChargedParticles_SourceObjects(locSourceChargedTracks);

	for(size_t loc_i = 0; loc_i < locSourceChargedTracks.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locSourceChargedTracks[loc_i]->dChargedTrackHypotheses[0];

		vector<const DTrackWireBased*>::iterator locIterator;
		for(locIterator = locUnusedWireBasedTracks.begin(); locIterator != locUnusedWireBasedTracks.end();)
		{
			if(locChargedTrackHypothesis->candidateid != (*locIterator)->candidateid)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedWireBasedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedTrackCandidates(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackCandidate*>& locUnusedTrackCandidates) const
{
	locUnusedTrackCandidates.clear();
	locEventLoop->Get(locUnusedTrackCandidates);

	deque<const DChargedTrack*> locSourceChargedTracks;
	locParticleCombo->Get_DetectedFinalChargedParticles_SourceObjects(locSourceChargedTracks);

	set<unsigned int> locUsedCandidateIndices;
	for(size_t loc_i = 0; loc_i < locSourceChargedTracks.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locSourceChargedTracks[loc_i]->dChargedTrackHypotheses[0];
		locUsedCandidateIndices.insert(locChargedTrackHypothesis->candidateid - 1); //id = index + 1
	}

	for(int loc_i = locUnusedTrackCandidates.size() - 1; loc_i >= 0; --loc_i)
	{
		if(locUsedCandidateIndices.find(loc_i) != locUsedCandidateIndices.end())
			locUnusedTrackCandidates.erase(locUnusedTrackCandidates.begin() + loc_i);
	}
}

void DAnalysisUtilities::Get_UnusedNeutralShowers(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralShower*>& locUnusedNeutralShowers) const
{
	locUnusedNeutralShowers.clear();
	locEventLoop->Get(locUnusedNeutralShowers, dShowerSelectionTag.c_str());

	deque<const DNeutralShower*> locSourceNeutralShowers;
	locParticleCombo->Get_DetectedFinalNeutralParticles_SourceObjects(locSourceNeutralShowers);

	for(size_t loc_i = 0; loc_i < locSourceNeutralShowers.size(); ++loc_i)
	{
		vector<const DNeutralShower*>::iterator locIterator;
		for(locIterator = locUnusedNeutralShowers.begin(); locIterator != locUnusedNeutralShowers.end();)
		{
			if(locSourceNeutralShowers[loc_i] != *locIterator)
				++locIterator;
			else
				locIterator = locUnusedNeutralShowers.erase(locIterator);
		}
	}
}

void DAnalysisUtilities::Get_UnusedNeutralParticles(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralParticle*>& locUnusedNeutralParticles) const
{
	locUnusedNeutralParticles.clear();
	locEventLoop->Get(locUnusedNeutralParticles, dShowerSelectionTag.c_str());

	deque<const DNeutralShower*> locSourceNeutralShowers;
	locParticleCombo->Get_DetectedFinalNeutralParticles_SourceObjects(locSourceNeutralShowers);

	for(size_t loc_i = 0; loc_i < locSourceNeutralShowers.size(); ++loc_i)
	{
		vector<const DNeutralParticle*>::iterator locIterator;
		for(locIterator = locUnusedNeutralParticles.begin(); locIterator != locUnusedNeutralParticles.end();)
		{
			if(locSourceNeutralShowers[loc_i] != (*locIterator)->dNeutralShower)
				++locIterator;
			else
				locIterator = locUnusedNeutralParticles.erase(locIterator);
		}
	}
}

void DAnalysisUtilities::Get_ThrownParticleSteps(JEventLoop* locEventLoop, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps) const
{
 	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	map<size_t, const DMCThrown*> locIDMap; //size_t is the myid
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		locIDMap[locMCThrowns[loc_i]->myid] = locMCThrowns[loc_i];

	locThrownSteps.clear();
	if(locMCThrowns.empty())
		return;
	locThrownSteps.push_back(pair<const DMCThrown*, deque<const DMCThrown*> >(NULL, deque<const DMCThrown*>() ) );

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(IsResonance(locMCThrowns[loc_i]->PID()))
			continue; //don't include resonances in DReaction!!

		if(locMCThrowns[loc_i]->PID() == Unknown)
			continue; //could be some weird pythia "resonance" like a diquark: just ignore them all

		//initial checks of parent id
		int locParentID = locMCThrowns[loc_i]->parentid;
		if(locParentID == 0) //photoproduced
		{
			locThrownSteps[0].second.push_back(locMCThrowns[loc_i]);
			continue;
		}
		if(locIDMap.find(locParentID) == locIDMap.end()) //produced from a particle that was not saved: spurious, don't save (e.g. product of BCAL shower)
			continue;

		//initial checks of parent pid
		Particle_t locParentPID = locIDMap[locParentID]->PID();
		bool locDoneFlag = false;
		while(((locParentPID == Unknown) || IsResonance(locParentPID)) && (!locDoneFlag))
		{
			//intermediate particle, continue towards the source
			locParentID = locIDMap[locParentID]->parentid; //parent's parent
			if(locParentID == 0) //photoproduced
			{
				locThrownSteps[0].second.push_back(locMCThrowns[loc_i]);
				locDoneFlag = true;
			}
			else if(locIDMap.find(locParentID) == locIDMap.end()) //produced from a particle that was not saved: spurious, don't save (e.g. product of BCAL shower)
				locDoneFlag = true;
			else
				locParentPID = locIDMap[locParentID]->PID();
		}

		if(Is_FinalStateParticle(locParentPID) == 1)
			continue; //e.g. parent is a final state particle (e.g. this is a neutrino from a pion decay)

		//check to see if the parent is already listed as a decaying particle //if so, add it to that step
		bool locListedAsDecayingFlag = false;
		for(size_t loc_j = 1; loc_j < locThrownSteps.size(); ++loc_j)
		{
			if(locThrownSteps[loc_j].first->myid != locParentID)
				continue;
			locThrownSteps[loc_j].second.push_back(locMCThrowns[loc_i]);
			locListedAsDecayingFlag = true;
			break;
		}
		if(locListedAsDecayingFlag)
			continue;

		//would add a new decay step, but first make sure that its parent is a decay product of a previous step
			//if the parent was not saved as a product, it may have been a decay product of a final state particle: don't save
		const DMCThrown* locThrownParent = locIDMap[locParentID];
		bool locFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locThrownSteps.size(); ++loc_j)
		{
			for(size_t loc_k = 0; loc_k < locThrownSteps[loc_j].second.size(); ++loc_k)
			{
				if(locThrownSteps[loc_j].second[loc_k] != locThrownParent)
					continue;
				locFoundFlag = true;
				break;
			}
			if(locFoundFlag)
				break;
		}
		if(!locFoundFlag)
			continue;

		//else add a new decay step and add this particle to it
		locThrownSteps.push_back(pair<const DMCThrown*, deque<const DMCThrown*> >(locThrownParent, deque<const DMCThrown*>(1, locMCThrowns[loc_i]) ));
	}

/*
cout << "THROWN STEPS: " << endl;
for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
{
	cout << ((loc_i == 0) ? 0 : locThrownSteps[loc_i].first->myid) << ": ";
	for(size_t loc_j = 0; loc_j < locThrownSteps[loc_i].second.size(); ++loc_j)
		cout << locThrownSteps[loc_i].second[loc_j]->myid << ", ";
	cout << endl;
}
*/
}

bool DAnalysisUtilities::Are_ThrownPIDsSameAsDesired(JEventLoop* locEventLoop, const deque<Particle_t>& locDesiredPIDs, Particle_t locMissingPID) const
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns, "FinalState");
	deque<Particle_t> locDesiredPIDs_Copy = locDesiredPIDs;

	bool locMissingPIDMatchedFlag = false;
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		Particle_t locPID = (Particle_t)(locMCThrowns[loc_i]->type);

		if((!locMissingPIDMatchedFlag) && (locMissingPID == locPID))
		{
			//matched missing
			locMissingPIDMatchedFlag = true;
			continue;
		}

		bool locPIDFoundFlag = false;
		for(deque<Particle_t>::iterator locIterator = locDesiredPIDs_Copy.begin(); locIterator != locDesiredPIDs_Copy.end(); ++locIterator)
		{
			if(*locIterator != locPID)
				continue;
			locDesiredPIDs_Copy.erase(locIterator);
			locPIDFoundFlag = true;
			break;
		}
		if(!locPIDFoundFlag)
			return false;
	}

	return (locDesiredPIDs_Copy.empty());
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DParticleCombo* locParticleCombo, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, Particle_t> > locSourceObjects;
	return Calc_MissingP4(locParticleCombo, 0, -1, deque<Particle_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DParticleCombo* locParticleCombo, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	return Calc_MissingP4(locParticleCombo, 0, -1, deque<Particle_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, Particle_t> > locSourceObjects;
	return Calc_MissingP4(locParticleCombo, locStepIndex, locUpToStepIndex, locUpThroughPIDs, locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	//NOTE: this routine assumes that the p4 of a charged decaying particle with a detached vertex is the same at both vertices!
	//assumes missing particle is not the beam particle
	if(locUseKinFitDataFlag && (locParticleCombo->Get_KinFitResults() == NULL))
		return Calc_MissingP4(locParticleCombo, locStepIndex, locUpToStepIndex, locUpThroughPIDs, locSourceObjects, false); //kinematic fit failed

	DLorentzVector locMissingP4;
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);

	const DKinematicData* locKinematicData = NULL;
	if(locStepIndex == 0)
	{
		//initial particle
		locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		locSourceObjects.insert(pair<const JObject*, Particle_t>(locKinematicData, locKinematicData->PID())); //want to use source objects for comparing
		if(locUseKinFitDataFlag) //kinfit
			locKinematicData = locParticleComboStep->Get_InitialParticle();
		locMissingP4 += locKinematicData->lorentzMomentum();

		//target particle
		locKinematicData = locParticleComboStep->Get_TargetParticle();
		if(locKinematicData != NULL)
		{
			locSourceObjects.insert(pair<const JObject*, Particle_t>(locKinematicData, locKinematicData->PID()));
			locMissingP4 += locKinematicData->lorentzMomentum();
		}
	}

	deque<const DKinematicData*> locParticles;
	if(!locUseKinFitDataFlag) //measured
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	else //kinfit
		locParticleComboStep->Get_FinalParticles(locParticles);

	for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
	{
		//DecayStepIndex: one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		if((locDecayStepIndex == -1) || (locDecayStepIndex == -3))
			continue; //missing particle or no blueprint

		Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
		if(int(locStepIndex) == locUpToStepIndex)
		{
			bool locPIDFoundFlag = false;
			for(deque<Particle_t>::iterator locIterator = locUpThroughPIDs.begin(); locIterator != locUpThroughPIDs.end(); ++locIterator)
			{
				if((*locIterator) != locPID)
					continue;
				locUpThroughPIDs.erase(locIterator);
				locPIDFoundFlag = true;
				break;
			}
			if(!locPIDFoundFlag)
				continue; //skip it: don't want to include it
		}

		if(locDecayStepIndex == -2) //detected
		{
			locMissingP4 -= locParticles[loc_j]->lorentzMomentum();
			locSourceObjects.insert(pair<const JObject*, Particle_t>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_j), locPID));
		}
		else //decaying-particle
			locMissingP4 += Calc_MissingP4(locParticleCombo, locDecayStepIndex, locUpToStepIndex, locUpThroughPIDs, locSourceObjects, locUseKinFitDataFlag); //p4 returned is already < 0
	}

	return locMissingP4;
}

DMatrixDSym DAnalysisUtilities::Calc_MissingP3Covariance(const DParticleCombo* locParticleCombo) const
{
	//uses measured data!
	return Calc_MissingP3Covariance(locParticleCombo, 0, -1, deque<Particle_t>());
}

DMatrixDSym DAnalysisUtilities::Calc_MissingP3Covariance(const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs) const
{
	//uses measured data!

	//NOTE: this routine assumes that the p4 of a charged decaying particle with a detached vertex is the same at both vertices!
	//assumes missing particle is not the beam particle

	//missing covariance is just sum of covariance matrices of all particles used in the calculation
		//because errors are uncorrelated: this doesn't work on kinfit data: just use kinfit matrix from missing particle then
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	DMatrixDSym locMissingCovarianceMatrix(3);
	locMissingCovarianceMatrix.Zero();

	const DKinematicData* locKinematicData = NULL;
	if(locStepIndex == 0)
	{
		//initial particle
		locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		DMatrixDSym locParticleCovarianceMatrix = locKinematicData->errorMatrix();
		locParticleCovarianceMatrix.ResizeTo(3, 3);
		locMissingCovarianceMatrix += locParticleCovarianceMatrix;
	}

	deque<const DKinematicData*> locParticles;
	locParticleComboStep->Get_FinalParticles_Measured(locParticles);

	for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
	{
		//DecayStepIndex: one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		if((locDecayStepIndex == -1) || (locDecayStepIndex == -3))
			continue; //missing particle or no blueprint

		Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
		if(int(locStepIndex) == locUpToStepIndex)
		{
			bool locPIDFoundFlag = false;
			for(deque<Particle_t>::iterator locIterator = locUpThroughPIDs.begin(); locIterator != locUpThroughPIDs.end(); ++locIterator)
			{
				if((*locIterator) != locPID)
					continue;
				locUpThroughPIDs.erase(locIterator);
				locPIDFoundFlag = true;
				break;
			}
			if(!locPIDFoundFlag)
				continue; //skip it: don't want to include it
		}

		if(locDecayStepIndex == -2) //detected
		{
			DMatrixDSym locParticleCovarianceMatrix = locParticles[loc_j]->errorMatrix();
			locParticleCovarianceMatrix.ResizeTo(3, 3);
			locMissingCovarianceMatrix += locParticleCovarianceMatrix;
		}
		else //decaying-particle
			locMissingCovarianceMatrix += Calc_MissingP3Covariance(locParticleCombo, locDecayStepIndex, locUpToStepIndex, locUpThroughPIDs);
	}

	return locMissingCovarianceMatrix;
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, Particle_t> > locSourceObjects;
	return Calc_FinalStateP4(locParticleCombo, locStepIndex, locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DParticleCombo* locParticleCombo, size_t locStepIndex, set<pair<const JObject*, Particle_t> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	if(locUseKinFitDataFlag && (locParticleCombo->Get_KinFitResults() == NULL))
		return Calc_FinalStateP4(locParticleCombo, locStepIndex, locSourceObjects, false); //kinematic fit failed

	DLorentzVector locFinalStateP4;
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	if(locParticleComboStep == NULL)
		return (DLorentzVector());

	deque<const DKinematicData*> locParticles;
	if(!locUseKinFitDataFlag) //measured
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	else //kinfit
		locParticleComboStep->Get_FinalParticles(locParticles);

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(locParticleComboStep->Is_FinalParticleMissing(loc_i))
			return (DLorentzVector()); //bad!
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
		{
			//measured results, or not constrained by kinfit (either non-fixed mass or excluded from kinfit)
			if((!locUseKinFitDataFlag) || (!IsFixedMass(locParticleComboStep->Get_FinalParticleID(loc_i))))
				locFinalStateP4 += Calc_FinalStateP4(locParticleCombo, locParticleComboStep->Get_DecayStepIndex(loc_i), locSourceObjects, locUseKinFitDataFlag);
			else //want kinfit results, and decaying particle p4 is constrained by kinfit
			{
				locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
				//still need source objects of decay products! dive down anyway, but ignore p4 result
				Calc_FinalStateP4(locParticleCombo, locParticleComboStep->Get_DecayStepIndex(loc_i), locSourceObjects, locUseKinFitDataFlag);
			}
		}
		else
		{
			locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
			locSourceObjects.insert(pair<const JObject*, Particle_t>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i), locParticles[loc_i]->PID()));
		}
	}
	return locFinalStateP4;
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinematicData, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir = (1.0/locKinematicData->momentum().Mag())*locKinematicData->momentum();
	DVector3 locPosition = locKinematicData->position();
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinFitParticle, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir(locKinFitParticle->Get_Momentum().Unit().X(),locKinFitParticle->Get_Momentum().Unit().Y(),locKinFitParticle->Get_Momentum().Unit().Z());
	DVector3 locPosition(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z());
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3& locDOCAVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	double locDOCA = Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
	locDOCAVertex = 0.5*(locPOCA1 + locPOCA2);
	return locDOCA;
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinFitParticle1, locKinFitParticle2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinematicData1, locKinematicData2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
  //originated from code by Jörn Langheinrich
  //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
  double locUnitDot = locUnitDir1.Dot(locUnitDir2);
  double locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
  double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point

  if(fabs(locDenominator) < 1.0e-15) //parallel
    locDistVertToInterDOCA1 = (locVertex2 - locVertex1).Dot(locUnitDir2)/locUnitDot; //the opposite
  else{
    double locA = (locVertex1 - locVertex2).Dot(locUnitDir1);
    double locB = (locVertex1 - locVertex2).Dot(locUnitDir2);
    locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
    locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
  }

  locPOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
  locPOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
  return (locPOCA1 - locPOCA2).Mag();
}

double DAnalysisUtilities::Calc_CrudeTime(const deque<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	DVector3 locMomentum;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
		locDeltaVertex = locPOCA - locParticles[loc_i]->position();
		locMomentum = locParticles[loc_i]->momentum();
		double locTime = locParticles[loc_i]->time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->energy()/(29.9792458*locMomentum.Mag2());
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

double DAnalysisUtilities::Calc_CrudeTime(const deque<const DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		double locTime = 0.0;
		double locE = locParticles[loc_i]->Get_ShowerEnergy();
		if((locParticles[loc_i]->Get_Charge() == 0) && (locE > 0.0))
		{
			double locMass = locParticles[loc_i]->Get_Mass();
			double locPMag = sqrt(locE*locE - locMass*locMass);
			TVector3 locPosition = locParticles[loc_i]->Get_Position();
			DVector3 locDPosition(locPosition.X(), locPosition.Y(), locPosition.Z());
			DVector3 locDeltaVertex = locDPosition - locCommonVertex;
			locTime = locParticles[loc_i]->Get_Time() + locDeltaVertex.Mag()*locE/(29.9792458*locPMag);
		}
		else
		{
			Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
			locDeltaVertex = locPOCA - DVector3(locParticles[loc_i]->Get_Position().X(),locParticles[loc_i]->Get_Position().Y(),locParticles[loc_i]->Get_Position().Z());
			DVector3 locMomentum(locParticles[loc_i]->Get_Momentum().X(),locParticles[loc_i]->Get_Momentum().Y(),locParticles[loc_i]->Get_Momentum().Z());
			locTime = locParticles[loc_i]->Get_Time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->Get_Energy()/(29.9792458*locMomentum.Mag2());
		}
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DTrackTimeBased*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DChargedTrackHypothesis*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DKinematicData*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const deque<const DKinFitParticle*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;

	if(locParticles.size() == 1)
	  return DVector3(locParticles[0]->Get_Position().X(), locParticles[0]->Get_Position().Y(), locParticles[0]->Get_Position().Z());

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

