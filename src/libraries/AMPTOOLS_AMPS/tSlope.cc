
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
	assert( args.size() == 1 );
	
	slope  = AmpParameter( args[0] );
  
	// need to register any free parameters so the framework knows about them
	registerParameter( slope );
}


complex< GDouble >
tSlope::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	
	TLorentzVector resonance = p1 + p2;
	TLorentzVector target(0, 0, 0, 0.938);

        //GDouble s = (beam + target).M2();
        GDouble t = (target - recoil).M2();

	GDouble W = exp(slope * t);
	
	return complex< GDouble > ( sqrt(fabs(W)) );
}

