
#include <iostream>
#include <stdlib.h>

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "particleType.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

const double ProductionMechanism::kPi = 3.14159;

using namespace std;

ProductionMechanism::ProductionMechanism( Recoil recoil, Type type, double slope ) :
m_type( type ),
m_lowMass( 0 ),
m_highMass( 0 ),
m_slope( slope ),
m_lastWeight( 1. )
{	
  kMproton=ParticleMass(Proton);
  kMneutron=ParticleMass(Neutron);

  switch( recoil ){
    // I'm sure the distinction between these doesn't matter!  
  case kProton:  m_recMass = kMproton; break; //old value: 0.9382
  case kNeutron: m_recMass = kMneutron; break; //old value: 0.9395
  default:       m_recMass = kMproton; break; //old value: 0.9382
  }
}

void
ProductionMechanism::setMassRange( double low, double high ){
	
	m_lowMass = low;
	m_highMass = high;
}

void 
ProductionMechanism::setGeneratorType( Type type ){
  
  m_type = type;
}


TLorentzVector
ProductionMechanism::produceResonance( const TLorentzVector& beam ){
	
	TLorentzVector target( 0, 0, 0, kMproton );
	
	TLorentzRotation lab2cmBoost( -( target + beam ).BoostVector() );
	TLorentzRotation cm2labBoost( ( target + beam ).BoostVector() );
	
	double cmEnergy = ( lab2cmBoost * ( target + beam ) ).E();
	double beamMomCM = cmMomentum( cmEnergy, beam.M(), target.M() );
  
  double t, tMax, resMass, resMomCM;

  do {
    resMass = generateMass();
    resMomCM  = cmMomentum( cmEnergy, resMass, m_recMass );
  
    tMax = 4. * beamMomCM * resMomCM;
    t = random( 0, tMax ); 
  } 
  while( random( 0., 1. ) > exp(-m_slope*t) );
	
	TVector3 resonanceMomCM;
	resonanceMomCM.SetMagThetaPhi( resMomCM,
                              acos( 1. - 2.*t/tMax ),
                              random( -kPi, kPi ) );
	
	TLorentzVector resonanceCM( resonanceMomCM, 
                               sqrt( resonanceMomCM.Mag2() +
                                    resMass * resMass ) );
	
	return cm2labBoost * resonanceCM;
}

void 
ProductionMechanism::addResonance( double mass, double width, double crossSec ){
  
  m_decGen.addChannel( m_bwGen.size(), crossSec );
  m_bwGen.push_back( BreitWignerGenerator( mass, width ) );
}

double
ProductionMechanism::generateMass(){
  
  if( m_type == kFlat ) return random( m_lowMass, m_highMass );
  
  double mass = 0;
  while( mass < m_lowMass || mass > m_highMass ){
    
    unsigned int channel = m_decGen();
    pair< double, double > bw = m_bwGen[channel]();
    
    mass = bw.first;
  }
  
  double prob = 0;
  for( unsigned int i = 0; i < m_bwGen.size(); ++i ){
    
    prob += m_bwGen[i].pdf( mass * mass ) * m_decGen.getProb( i );
  }
  
  // put in the factor of mass so resulting weights can be applied to 
  // obtain a distribution that is flat in mass instead of flat in s
  // (to get weights for reweighting flat in s remove the mass)
  
  m_lastWeight = 1 / ( prob * mass );
  
  return mass;
}

double
ProductionMechanism::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
ProductionMechanism::random( double low, double hi ) const {
	
	return( ( hi - low ) * drand48() + low );
}


