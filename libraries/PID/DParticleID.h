// $Id$
//
//    File: DParticleID.h
// Created: Mon Feb 28 13:47:49 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_
#define _DParticleID_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <DVector3.h>
#include "HDGEOMETRY/DRootGeom.h"
#include <TRACKING/DTrackTimeBased_factory.h>
class DTrackTimeBased;
class DCDCTrackHit;

class DParticleID:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DParticleID);
 
  // Constructor and destructor
  DParticleID(JEventLoop *loop); // require JEventLoop in constructor
  virtual ~DParticleID();

  class dedx_t{
  public:
    dedx_t(double dE,double dx, double p):dE(dE),dx(dx),p(p){}
      double dE; // energy loss in layer
      double dx; // path length in layer
      double p;  // momentum at this dE/dx measurement
      
  };


  virtual jerror_t GetdEdxChiSq(const DTrackTimeBased *track,double &dEdx,
				unsigned int &num,double &chi2)=0;
  double GetdEdxSigma(double num_hits,double p,double mass,
		      double mean_path_length);
  double GetMostProbabledEdx(double p,double mass,double dx);
  jerror_t GetdEdx(const DTrackTimeBased *track,vector<dedx_t>&dEdx_list);
  jerror_t CalcdEdxHit(const DVector3 &mom,
		       const DVector3 &pos,
		       const DCDCTrackHit *hit,
		       pair <double,double> &dedx);
  jerror_t GroupTracks(vector<const DTrackTimeBased *> &tracks,
		       vector<vector<const DTrackTimeBased*> >&grouped_tracks);



 private: 
  //< DGeometry pointer used to access materials through calibDB maps for eloss
  const DRootGeom *RootGeom;                                 
 
  
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID();
  

  // gas material properties
  double mKRhoZoverAGas,mRhoZoverAGas,mLnIGas;

		
};

#endif // _DParticleID_

