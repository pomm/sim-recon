

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPiSoding.h"

TwoPiSoding::TwoPiSoding( const vector< string >& args ) :
UserAmplitude< TwoPiSoding >( args )
{
  
  assert( args.size() == 16);
	
  m_mass0 = AmpParameter( args[0] );
  m_width0 = AmpParameter( args[1] );
  m_orbitL = atoi( args[2].c_str() );
  m_daughters = pair< string, string >( args[3], args[4] );
  m_soding = AmpParameter( args[5] );
  m_slope = AmpParameter( args[6] );

  rho000  = AmpParameter( args[7] );
  rho100  = AmpParameter( args[8] );
  rho1m10 = AmpParameter( args[9] );
  
  rho111  = AmpParameter( args[10] );
  rho001  = AmpParameter( args[11] );
  rho101  = AmpParameter( args[12] );
  rho1m11 = AmpParameter( args[13] );
  
  rho102  = AmpParameter( args[14] );
  rho1m12 = AmpParameter( args[15] );

  // need to register any free parameters so the framework knows about them
  registerParameter( m_mass0 );
  registerParameter( m_width0 );
  registerParameter( m_soding );
  registerParameter( m_slope );
  
  registerParameter( rho000 );
  registerParameter( rho100 );
  registerParameter( rho1m10 );
  
  registerParameter( rho111 );
  registerParameter( rho001 );
  registerParameter( rho101 );
  registerParameter( rho1m11 );
  
  registerParameter( rho102 );
  registerParameter( rho1m12 );

  // make sure the input variables look reasonable
  assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );
}

complex< GDouble >
TwoPiSoding::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp;
  
  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
    
    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P1 += Ptemp;
    Ptot += Ptemp;
  }
  
  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
    
    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P2 += Ptemp;
    Ptot += Ptemp;
  }
  
  GDouble mass  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();
  
  // assert positive breakup momenta     
  GDouble q0 = fabs( breakupMomentum(m_mass0, mass1, mass2) );
  GDouble q  = fabs( breakupMomentum(mass,    mass1, mass2) );
  
  GDouble F0 = barrierFactor(q0, m_orbitL);
  GDouble F  = barrierFactor(q,  m_orbitL);
  
  GDouble width = m_width0*(m_mass0/mass)*(q/q0)*((F*F)/(F0*F0));
  
  // this first factor just gets normalization right for BW's that have
  // no additional s-dependence from orbital L
  complex<GDouble> bwtop( sqrt( m_mass0 * m_width0 / 3.1416 ), 0.0 );
  complex<GDouble> bwbottom( ( m_mass0*m_mass0 - mass*mass ) ,
                           -1.0 * ( m_mass0 * width ) );
  
  complex<GDouble> bw = F * bwtop / bwbottom;

  // compute Soding factor
  complex<GDouble> sodingtop( m_mass0*m_mass0 - mass*mass );
  complex<GDouble> sodingbottom( m_mass0*m_mass0 - mass*mass, -1.0 * m_mass0*m_width0 );
  
  complex<GDouble> soding = sodingtop / sodingbottom;

  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  
  TLorentzVector resonance = p1 + p2;
  GDouble t = (beam - resonance).M2();

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
  
  // vector meson production from K. Schilling et. al.
  GDouble Pgamma = 0.0;
  GDouble Egamma = beam.E();
  if(Egamma < 3.0)
    Pgamma = 0.6 * TMath::Gaus(Egamma, 3.0, 0.3);
  
  GDouble W = 0.5*(1. - rho000) + 0.5*(3.*rho000 - 1.)*cosTheta*cosTheta - sqrt(2.)*rho100*sin2Theta*cos(phi) - rho1m10*sinSqTheta*cos(2.*phi);
  
  W -= Pgamma*cos(2.*Phi) * (rho111*sinSqTheta + rho001*cosTheta*cosTheta - sqrt(2.)*rho101*sin2Theta*cos(phi) - rho1m11*sinSqTheta*cos(2.*phi));
  
  W -= Pgamma*sin(2.*Phi) * (sqrt(2.)*rho102*sin2Theta*sin(phi) + rho1m12*sinSqTheta*sin(2.*phi));
  
  W *= 3./(4.*PI);
  
  GDouble tot = norm(bw) * W * exp( m_slope * t ); 
  tot += 2.0 * m_soding * real(conj(soding)*bw) * exp( m_slope * t/2.0 );
  tot += m_soding*m_soding * norm(soding);

  return( sqrt( fabs(tot) ) );
}

void
TwoPiSoding::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

#ifdef GPU_ACCELERATION
void
TwoPiSoding::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  // use integers to endcode the string of daughters -- one index in each
  // decimal place
  
  int daught1 = atoi( m_daughters.first.c_str() );
  int daught2 = atoi( m_daughters.second.c_str() );
  
  GPUTwoPiSoding_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, 
                       m_mass0, m_width0, m_orbitL, daught1, daught2 );

}
#endif //GPU_ACCELERATION

