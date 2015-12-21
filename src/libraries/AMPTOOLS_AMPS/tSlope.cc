
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/tSlope.h"

tSlope::tSlope( const vector< string >& args ) :
UserAmplitude< tSlope >( args )
{
	assert( args.size() > 1 );
	
	slope  = AmpParameter( args[0] );

	// need to register any free parameters so the framework knows about them
	registerParameter( slope );

	// set children of resonance to sum
	for(unsigned int i = 1; i < args.size(); i++) {
	  resonanceChildren.push_back( atoi( args[i].data() ) );
	}
}


complex< GDouble >
tSlope::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 

	// sum p4 from children of resonace
	TLorentzVector resonance;
	for(unsigned int i = 0; i < resonanceChildren.size(); i++) {
	  TLorentzVector particle ( pKin[resonanceChildren[i]][1], pKin[resonanceChildren[i]][2], pKin[resonanceChildren[i]][3], pKin[resonanceChildren[i]][0] ); 
	  resonance += particle;
	}

        GDouble t = (beam - resonance).M2();

	GDouble W = exp(slope * t);
	
	return complex< GDouble > ( sqrt(fabs(W)) );
}

