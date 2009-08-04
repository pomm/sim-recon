// $Id$
//
//    File: DReferenceTrajectory.h
// Created: Wed Jul 19 13:42:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DReferenceTrajectory_
#define _DReferenceTrajectory_

#include <vector>
using std::vector;

#include <HDGEOMETRY/DGeometry.h>
#include <DVector3.h>

#include <JANA/jerror.h>

#include <DCoordinateSystem.h>

class DMagneticFieldMap;
class DTrackCandidate;
class DRootGeom;

class DReferenceTrajectory{
	
	/// This class is a utility class used by the TRACKING package. It
	/// is used to swim a particle through the (inhomogeneous) magnetic
	/// field, accounting for energy loss in the geometry, and recording
	/// each step so that derivatives and distances can be calculated.
	/// Because this uses the coordinates defined in the reference 
	/// trajectory method (B. Mecking Nucl. Instr. and Methods. 
	/// 203 (1982) 299-305), the angles needed to rotate into the
	/// lab frame are saved as well. This used the DMagneticFieldStepper
	/// class for swimming through the field.

	public:

		enum direction_t{
			kForward,
			kBackward
		};

		class swim_step_t:public DCoordinateSystem{
			public:
				DVector3 mom;
				double Ro;
				double s; // distance along RT
				double dP;
				
				// The following are used to calculate the covariance matrix for MULS
				double itheta02;		// running sum of MULS angle theta_0 squared
				double itheta02s;		// ditto but times s
				double itheta02s2;	// ditto but times s^2
		};

		DReferenceTrajectory(const DMagneticFieldMap *
									, double q=1.0
									, swim_step_t *swim_steps=NULL
									, int max_swim_steps=0
									, double step_size=-1.0);
		DReferenceTrajectory(const DReferenceTrajectory& rt);
		DReferenceTrajectory& operator=(const DReferenceTrajectory& rt);
		virtual ~DReferenceTrajectory();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DReferenceTrajectory";}
		
		double DistToRT(double x, double y, double z){return DistToRT(DVector3(x,y,z));}
		double DistToRT(DVector3 hit);
		double DistToRT(const DCoordinateSystem *wire, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, double *s=NULL);
		double DistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		double DistToRTBruteForce(const DCoordinateSystem *wire, const swim_step_t *step, double *s=NULL);
		double Straw_dx(const DCoordinateSystem *wire, double radius);
		swim_step_t* FindClosestSwimStep(const DCoordinateSystem *wire, int *istep_ptr=NULL);
		void Swim(const DVector3 &pos, const DVector3 &mom, double q=-1000.0, double smax=2000.0, const DCoordinateSystem *wire=NULL);
		int InsertSteps(const swim_step_t *start_step, double delta_s, double step_size=0.02); 
		DVector3 GetLastDOCAPoint(void);
		void GetLastDOCAPoint(DVector3 &pos, DVector3 &mom);
		double GetLastDistAlongWire(void){return last_dist_along_wire;}
		void SetStepSize(double step_size){this->step_size=step_size;}
		void SetDRootGeom(const DRootGeom *RootGeom){this->RootGeom = RootGeom;}
		void SetDGeometry(const DGeometry *geom){this->geom = geom;}
		const DRootGeom* GetDRootGeom(void){return RootGeom;}
		const DGeometry* GetDGeometry(void){return geom;}
		double GetMass(void) const {return mass;}
		void SetMass(double mass){this->mass = mass;}
		void SetPLossDirection(direction_t direction){ploss_direction=direction;}
		direction_t GetPLossDirection(void){return ploss_direction;}
		inline double dPdx(double ptot, double A, double Z, double density);
		inline double dPdx(double ptot, double rhoZ_overA, double rhoZ_overA_logI);

		const swim_step_t *GetLastSwimStep(void){return last_swim_step;}
		swim_step_t *swim_steps;
		int Nswim_steps;
		float q;

	protected:
	
		int max_swim_steps;
		bool own_swim_steps;
		int dist_to_rt_depth;
		double step_size;
		const DMagneticFieldMap *bfield;
		const DRootGeom *RootGeom;
		const DGeometry *geom;
		direction_t ploss_direction;
		
		double last_phi;							///< last phi found in DistToRT
		const swim_step_t* last_swim_step;	///< last swim step used in DistToRT
		double last_dist_along_wire;
		double last_dz_dphi;
		
		double mass;
	
	private:
		DReferenceTrajectory(){} // force use of constructor with arguments.

};

#endif // _DReferenceTrajectory_

