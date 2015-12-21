
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiRhoJPAC.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

// FORTRAN routines
extern "C"{
double xnorm_(float* s, float* t, float* xm, float* xnorm);
dcomplex rho_h_(int* k, int* l, int* lp, float* s, float* t, float* xm);
extern struct{
	dcomplex param[4];
} pars_;
extern struct{
	double b[4];
} slopes_;
};

// Wrapper function for differential cross section
GDouble Norm(double s, double t, double xm)
{
        float xx = s;
	float yy = t;
	float zz = xm;
	float xnorm = 0.;
	xnorm_(&xx, &yy, &zz, &xnorm);

	if(xnorm!=0) {
		if(xnorm!=xnorm) {
			//cout<<endl<<"s="<<xx<<" t="<<yy<<" xm="<<zz<<endl;
			//cout<<xnorm<<endl;
			return 0.;
		}
	}

	return xnorm;
}

// Wrapper function for SDME
GDouble rhoH(int k, int l, int lp, double s, double t, double xm, bool isReal)
{
	int x = k;
	int y = l;
	int z = lp;
	
        float xx = s;
	float yy = t;
	float zz = xm;
	dcomplex rho_h = rho_h_(&x, &y, &z, &xx, &yy, &zz);
	
	double rho = 0;
	if(isReal) rho = rho_h.dr;
	else rho = rho_h.di;

	if(rho!=rho) {
		//cout<<endl<<"s="<<xx<<" t="<<yy<<" xm="<<zz<<endl;
		//cout<<rho<<endl;
		return 0.;
	}

	return rho;
}

TwoPiRhoJPAC::TwoPiRhoJPAC( const vector< string >& args ) :
UserAmplitude< TwoPiRhoJPAC >( args )
{
	assert( args.size() == 8 );
	
	param0  = AmpParameter( args[0] );
	param1  = AmpParameter( args[1] );
	param2  = AmpParameter( args[2] );
	param3  = AmpParameter( args[3] );

	b0  = AmpParameter( args[4] );
	b1  = AmpParameter( args[5] );
	b2  = AmpParameter( args[6] );	
	b3  = AmpParameter( args[7] );

	// need to register any free parameters so the framework knows about them
	registerParameter( param0 );
	registerParameter( param1 );
	registerParameter( param2 );
	registerParameter( param3 );

	registerParameter( b0 );
	registerParameter( b1 );
	registerParameter( b2 );
	registerParameter( b3 );
	
}


complex< GDouble >
TwoPiRhoJPAC::calcAmplitude( GDouble** pKin ) const {
  
	// set parameters for fortran code
	pars_.param[0].dr = param0;
	pars_.param[1].dr = param1;
	pars_.param[2].dr = param2;
	pars_.param[3].dr = param3;
	pars_.param[0].di = 0.;
	pars_.param[1].di = 0.;
	pars_.param[2].di = 0.;
	pars_.param[3].di = 0.;
	slopes_.b[0] = b0;
	slopes_.b[1] = b1;
	slopes_.b[2] = b2;
	slopes_.b[3] = b3;

	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
	
	TLorentzVector resonance = p1 + p2;
	TLorentzRotation resonanceBoost( -resonance.BoostVector() );
	TLorentzVector beam_res = resonanceBoost * beam;
	TLorentzVector recoil_res = resonanceBoost * recoil;
	TLorentzVector p1_res = resonanceBoost * p1;
	
	TVector3 z = -recoil_res.Vect().Unit();
	TVector3 y = beam_res.Vect().Cross(z).Unit();
	TVector3 x = y.Cross(z).Unit();
	TVector3 angles(   (p1_res.Vect()).Dot(x),
			   (p1_res.Vect()).Dot(y),
			   (p1_res.Vect()).Dot(z) );
	
	GDouble cosTheta = angles.CosTheta();
	GDouble sinSqTheta = sin(angles.Theta())*sin(angles.Theta());
	GDouble sin2Theta = sin(2.*angles.Theta());
	GDouble phi = angles.Phi();

	TVector3 eps(1.0, 0.0, 0.0); // beam linear polarization vector
	GDouble Phi = atan2(y.Dot(eps), beam_res.Vect().Unit().Dot(eps.Cross(y)));
	
	TLorentzVector target(0, 0, 0, 0.938);
	GDouble s = (beam + target).M2();
	GDouble t = (target - recoil).M2();
	GDouble xm = resonance.M();

	// vector meson production from K. Schilling et. al.
	GDouble Pgamma = 0.5;
	
	//rho(L,L',k,s,t,xm)
	GDouble N = Norm(s,t,xm);
	GDouble rho000 = rhoH(0,0,0,s,t,xm,1)/N;
	GDouble W = 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta*cosTheta - sqrt(2.)*rhoH(0,1,0,s,t,xm,1)/N*sin2Theta*cos(phi) - rhoH(0,1,-1,s,t,xm,1)/N*sinSqTheta*cos(2.*phi);
	W -= Pgamma*cos(2.*Phi) * (rhoH(1,1,1,s,t,xm,1)/N*sinSqTheta + rhoH(1,0,0,s,t,xm,1)/N*cosTheta*cosTheta - sqrt(2.)*rhoH(1,1,0,s,t,xm,1)/N*sin2Theta*cos(phi) - rhoH(1,1,-1,s,t,xm,1)/N*sinSqTheta*cos(2.*phi));
        W -= Pgamma*sin(2.*Phi) * (sqrt(2.)*rhoH(2,1,0,s,t,xm,0)/N*sin2Theta*sin(phi) + rhoH(2,1,-1,s,t,xm,0)/N*sinSqTheta*sin(2.*phi));
        W *= 3./(4.*PI);

	// amplitude from JPAC model
	W *= N;

	if(N==0) W=0;

	//if(xm > 0.7775 && xm < 0.7776 && s>6.34 && s<6.341) {
	//	cout<<s<<" "<<t<<" "<<xm<<endl;
	//	cout<<Norm(s,t,xm)<<endl;
	//	cout<<rhoH(2,1,-1,s,t,xm,0)/N<<" "<<rhoH(1,1,-1,s,t,xm,1)/N<<endl;
	//}

	//cout<<W<<endl<<endl;

	return complex< GDouble > ( sqrt(fabs(W)) );
}

