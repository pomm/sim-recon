/*
 *  GammaPToNPartP_Phasespace.cc
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "TLorentzVector.h"

#include "AMPTOOLS_MCGEN/GammaPToNPartP_Phasespace.h"

// FORTRAN routines
extern "C"{
void cobrems_(float* Emax, float* Epeak, float* emitmr, float* radt, float* Dist, float* collDiam, int* doPolFluxfloat);
float dntdx_(float* x);
float dnidx_(float* x);
};

// Wrapper function for total
double dNtdx(double x)
{
        float xx = x;
        return (double)dntdx_(&xx);
}

double dNidx(double x)
{
        float xx = x;
        return (double)dnidx_(&xx);
}

GammaPToNPartP_Phasespace::GammaPToNPartP_Phasespace( vector<double> &ChildMass, float beamMaxE, float beamPeakE, float beamLowE, float beamHighE) : 
m_target( 0, 0, 0, 0.938 ),
m_childMass( ChildMass ) {

  m_Npart = ChildMass.size();
  assert(m_Npart>1);
 
  // Initialize coherent brem table
  float Emax =  beamMaxE;
  float Epeak = beamPeakE;
  float Elow = beamLowE;
  float Ehigh = beamHighE;

  int doPolFlux=0;  // want total flux (1 for polarized flux)
  float emitmr=20.e-9; // electron beam emittance
  float radt=50.e-6; // radiator thickness in m
  float collDiam=0.005; // meters
  float Dist = 76.0; // meters
  cobrems_(&Emax, &Epeak, &emitmr, &radt, &Dist, &collDiam, &doPolFlux);

  // Create histogram
  cobrem_vs_E = new TH1D("cobrem_vs_E", "Coherent Bremstrahlung vs. E_{#gamma}", 1000, Elow, Ehigh);
  
  // Fill histogram
  for(int i=1; i<=cobrem_vs_E->GetNbinsX(); i++){
	  double x = cobrem_vs_E->GetBinCenter(i)/Emax;
	  double y = 0;
	  if(Epeak<Elow) y = dNidx(x);
	  else y = dNtdx(x);
	  cobrem_vs_E->SetBinContent(i, y);
  }

}

Kinematics* 
GammaPToNPartP_Phasespace::generate(){

  double beamE = cobrem_vs_E->GetRandom();
  m_beam.SetPxPyPzE(0,0,beamE,beamE);
  TLorentzVector cm = m_beam + m_target;

  Double_t masses[m_Npart];
  for(unsigned int i=0; i<m_Npart; i++)
	  masses[i] = m_childMass[i];

  TGenPhaseSpace phsp;
  phsp.SetDecay(cm,m_Npart,masses);

  double phsp_wt_max = phsp.GetWtMax();
  double genWeight;
  do {
     genWeight = phsp.Generate();
  }
  while( random(0., phsp_wt_max) >= genWeight || genWeight != genWeight);

  vector< TLorentzVector > allPart;
  allPart.push_back( m_beam );
  for(unsigned int i=0; i<m_Npart; i++) {
	  allPart.push_back( *phsp.GetDecay(i) );
  }
 
  return new Kinematics( allPart, genWeight );
}

double
GammaPToNPartP_Phasespace::random( double low, double hi ) const {

        return( ( hi - low ) * drand48() + low );
}

