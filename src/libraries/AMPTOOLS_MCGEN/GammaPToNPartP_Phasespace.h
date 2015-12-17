#if !defined(GAMMAPTONPARTP_PHASESPACE)
#define GAMMAPTONPARTP_PHASESPACE

/*
 *  GammaPToNPartP_Phasespace.h
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 1/22/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"

#include "IUAmpTools/Kinematics.h"

class Kinematics;
class AmpVecs;

class GammaPToNPartP_Phasespace {
  
public:
  
  GammaPToNPartP_Phasespace( vector<double> &ChildMass, float beamMaxE, float beamPeakE, float beamLowE, float beamHighE);
  
  Kinematics* generate();
  
private:
  
  TLorentzVector m_beam;
  TLorentzVector m_target;
  
  vector< double > m_childMass;
  unsigned int m_Npart;
  
  TH1D *cobrem_vs_E;
  
  double random( double low, double hi ) const;

};

#endif
