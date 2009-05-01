// $Id$
//
//    File: DReferenceTrajectory.cc
// Created: Wed Jul 19 13:42:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <memory>

#include <DVector3.h>
using namespace std;

#include "DReferenceTrajectory.h"
#include "DTrackCandidate.h"
#include "DMagneticFieldStepper.h"
#include "HDGEOMETRY/DRootGeom.h"

struct StepStruct {DReferenceTrajectory::swim_step_t steps[256];};

//---------------------------------
// DReferenceTrajectory    (Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DMagneticFieldMap *bfield
														, double q
														, swim_step_t *swim_steps
														, int max_swim_steps
														, double step_size)
{
	// It turns out that the greatest bottleneck in speed here comes from
	// allocating/deallocating the large block of memory required to hold
	// all of the trajectory info. The preferred way of calling this is 
	// with a pointer allocated once at program startup. This code block
	// though allows it to be allocated here if necessary.
	if(!swim_steps){
		own_swim_steps = true;
		this->max_swim_steps = 50000;
		this->swim_steps = new swim_step_t[this->max_swim_steps];
	}else{
		own_swim_steps = false;
		this->max_swim_steps = max_swim_steps;
		this->swim_steps = swim_steps;
	}

	// Copy some values into data members
	this->q = q;
	this->step_size = step_size;
	this->bfield = bfield;
	this->Nswim_steps = 0;
	this->dist_to_rt_depth = 0;

	this->RootGeom=NULL;
}

//---------------------------------
// DReferenceTrajectory    (Copy Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DReferenceTrajectory& rt)
{
	/// The copy constructor will always allocate its own memory for the
	/// swim steps and set its internal flag to indicate that is owns them
	/// regardless of the owner of the source trajectory's.

	this->Nswim_steps = rt.Nswim_steps;
	this->q = rt.q;
	this->max_swim_steps = rt.max_swim_steps;
	this->own_swim_steps = true;
	this->step_size = rt.step_size;
	this->bfield = rt.bfield;
	this->last_phi = rt.last_phi;
	this->last_swim_step = rt.last_swim_step;
	this->last_dist_along_wire = rt.last_dist_along_wire;
	this->last_dz_dphi = rt.last_dz_dphi;
	this->RootGeom = rt.RootGeom;
	this->dist_to_rt_depth = 0;

	this->swim_steps = new swim_step_t[this->max_swim_steps];
	for(int i=0; i<Nswim_steps; i++)swim_steps[i] = rt.swim_steps[i];
}

//---------------------------------
// operator=               (Assignment operator)
//---------------------------------
DReferenceTrajectory& DReferenceTrajectory::operator=(const DReferenceTrajectory& rt)
{
	/// The assignment operator will always make sure the memory allocated
	/// for the swim_steps is owned by the object being copied into.
	/// If it already owns memory of sufficient size, then it will be
	/// reused. If it owns memory that is too small, it will be freed and
	/// a new block allocated. If it does not own its swim_steps coming
	/// in, then it will allocate memory so that it does own it on the
	/// way out.
	
	if(&rt == this)return *this; // protect against self copies

	// Free memory if block is too small
	if(own_swim_steps==true && max_swim_steps<rt.Nswim_steps){
		delete[] swim_steps;
		swim_steps=NULL;
	}
	
	// Forget memory block if we don't currently own it
	if(!own_swim_steps){
		swim_steps=NULL;
	}

	this->Nswim_steps = rt.Nswim_steps;
	this->q = rt.q;
	this->max_swim_steps = rt.max_swim_steps;
	this->own_swim_steps = true;
	this->step_size = rt.step_size;
	this->bfield = rt.bfield;
	this->last_phi = rt.last_phi;
	this->last_swim_step = rt.last_swim_step;
	this->last_dist_along_wire = rt.last_dist_along_wire;
	this->last_dz_dphi = rt.last_dz_dphi;
	this->RootGeom = rt.RootGeom;
	this->dist_to_rt_depth = rt.dist_to_rt_depth;

	// Allocate memory if needed
	if(swim_steps==NULL)this->swim_steps = new swim_step_t[this->max_swim_steps];

	// Copy swim steps
	for(int i=0; i<Nswim_steps; i++)swim_steps[i] = rt.swim_steps[i];
	
	return *this;
}

//---------------------------------
// ~DReferenceTrajectory    (Destructor)
//---------------------------------
DReferenceTrajectory::~DReferenceTrajectory()
{
	if(own_swim_steps){
		delete[] swim_steps;
	}
}

//---------------------------------
// Swim
//---------------------------------
void DReferenceTrajectory::Swim(const DVector3 &pos, const DVector3 &mom, double q, double smax)
{
	/// (Re)Swim the trajectory starting from pos with momentum mom.
	/// This will use the charge and step size (if given) passed to
	/// the constructor when the object was created. It will also
	/// (re)use the sim_step buffer, replacing it's contents.
	if(fabs(q)>10)
		q = this->q;
	else
		this->q = q;

	DMagneticFieldStepper stepper(bfield, q, &pos, &mom);
	if(step_size>0.0)stepper.SetStepSize(step_size);

	// Step until we hit a boundary (don't track more than 20 meters)
	swim_step_t *swim_step = this->swim_steps;
	Nswim_steps = 0;
	double itheta02 = 0.0;
	double itheta02s = 0.0;
	double itheta02s2 = 0.0;
	swim_step_t *last_step=NULL;
	for(double s=0; fabs(s)<smax; Nswim_steps++, swim_step++){

		if(Nswim_steps>=this->max_swim_steps){
			cerr<<__FILE__<<":"<<__LINE__<<" Too many steps in trajectory. Truncating..."<<endl;
			break;
		}

		stepper.GetDirs(swim_step->sdir, swim_step->tdir, swim_step->udir);
		stepper.GetPosition(swim_step->origin);
		stepper.GetMomentum(swim_step->mom);
		swim_step->Ro = stepper.GetRo();
		swim_step->s = s;

		// Add material if RootGeom is not NULL
		if(RootGeom){
			double density, A, Z, X0;
			RootGeom->FindMat(swim_step->origin,density,A,Z,X0);
			if(X0>0.0){
				double delta_s = s;
				if(last_step)delta_s -= last_step->s;
				double radlen = delta_s/X0;
				if(radlen>1.0E-5){ // PDG 2008 pg 271, second to last paragraph
					double p = swim_step->mom.Mag();
					double beta = p/sqrt(p*p + 0.136*0.136); // assume pion mass
					double theta0 = 0.0136/(p*beta)*sqrt(radlen)*(1.0+0.038*log(radlen)); // From PDG 2008 eq 27.12
					double theta02 = theta0*theta0;
					itheta02 += theta02;
					itheta02s += s*theta02;
					itheta02s2 += s*s*theta02;
				}
			}
			last_step = swim_step;
		}
		swim_step->itheta02 = itheta02;
		swim_step->itheta02s = itheta02s;
		swim_step->itheta02s2 = itheta02s2;
		

		// Exit loop if we leave the tracking volume
		if(swim_step->origin.Perp()>65.0){break;} // ran into BCAL
		if(swim_step->origin.Z()>650.0){break;} // ran into FCAL
		if(swim_step->origin.Z()<-50.0){break;} // ran into UPV

		// Swim to next
		s += stepper.Step(NULL);
	}

	// OK. At this point the positions of the trajectory in the lab
	// frame have been recorded along with the momentum of the
	// particle and the directions of reference trajectory
	// coordinate system at each point.
}


//---------------------------------
// InsertSteps
//---------------------------------
int DReferenceTrajectory::InsertSteps(const swim_step_t *start_step, double delta_s, double step_size)
{
	/// Insert additional steps into the reference trajectory starting
	/// at start_step and swimming for at least delta_s by step_size
	/// sized steps. Both delta_s and step_size are in centimeters.
	/// If the value of delta_s is negative then the particle's momentum
	/// and charge are reversed before swimming. This could be a problem
	/// if energy loss is implemented.
	
	if(!start_step)return -1;

	// We do this by creating another, temporary DReferenceTrajectory object
	// on the stack and swimming it.
	DVector3 pos = start_step->origin;
	DVector3 mom = start_step->mom;
	double my_q = q;
	int direction = +1;
	if(delta_s<0.0){
		mom *= -1.0;
		my_q = -q;
		direction = -1;
	}
	
	// Here I allocate the steps using an auto_ptr so I don't have to mess with
	// deleting them at all of the possible exists. The problem with auto_ptr
	// is it can't handle arrays so it has to be wrapped in a struct.
	auto_ptr<StepStruct> steps_aptr(new StepStruct);
	DReferenceTrajectory::swim_step_t *steps = steps_aptr->steps;
	DReferenceTrajectory rt(bfield , my_q , steps , 256);
	rt.SetStepSize(step_size);
	rt.Swim(pos, mom, my_q, fabs(delta_s));
	if(rt.Nswim_steps==0)return 1;

	// Check that there is enough space to add these points
	if((Nswim_steps+rt.Nswim_steps)>max_swim_steps){
		//_DBG_<<"Not enough swim steps available to add new ones! Max="<<max_swim_steps<<" had="<<Nswim_steps<<" new="<<rt.Nswim_steps<<endl;
		return 2;
	}
	
	// At this point, we may have swum forward or backwards so the points
	// will need to be added either before start_step or after it. We also
	// may want to replace an old step that overlaps our high density steps
	// since they are presumably more accurate. Find the indexes of the
	// existing steps that the new steps will be inserted between.
	double sdiff = rt.swim_steps[rt.Nswim_steps-1].s;
	double s1 = start_step->s;
	double s2 = start_step->s + (double)direction*sdiff;
	double smin = s1<s2 ? s1:s2;
	double smax = s1<s2 ? s2:s1;
	int istep_start = 0;
	int istep_end = 0;
	for(int i=0; i<Nswim_steps; i++){
		if(swim_steps[i].s <  smin)istep_start = i;
		if(swim_steps[i].s <= smax)istep_end = i+1;
	}
	
	// istep_start and istep_end now point to the steps we want to keep.
	// All steps between them (exclusive) will be overwritten. Note that 
	// the original start_step should be in the "overwrite" range since 
	// it is included already in the new trajectory.
	int steps_to_overwrite = istep_end - istep_start - 1;
	int steps_to_shift = rt.Nswim_steps - steps_to_overwrite;
	
	// Shift the steps down (or is it up?) starting with istep_end.
	for(int i=Nswim_steps-1; i>=istep_end; i--)swim_steps[i+steps_to_shift] = swim_steps[i];
	
	// Copy the new steps into this object
	double s_0 = start_step->s;
	double itheta02_0 = start_step->itheta02;
	double itheta02s_0 = start_step->itheta02s;
	double itheta02s2_0 = start_step->itheta02s2;
	for(int i=0; i<rt.Nswim_steps; i++){
		int index = direction>0 ? (istep_start+1+i):(istep_start+1+rt.Nswim_steps-1-i);
		swim_steps[index] = rt.swim_steps[i];
		swim_steps[index].s = s_0 + (double)direction*swim_steps[index].s;
		swim_steps[index].itheta02 = itheta02_0 + (double)direction*swim_steps[index].itheta02;
		swim_steps[index].itheta02s = itheta02s_0 + (double)direction*swim_steps[index].itheta02s;
		swim_steps[index].itheta02s2 = itheta02s2_0 + (double)direction*swim_steps[index].itheta02s2;
		if(direction<0.0){
			swim_steps[index].sdir *= -1.0;
			swim_steps[index].tdir *= -1.0;
		}
	}
	Nswim_steps += rt.Nswim_steps-steps_to_overwrite;

	// Note that the above procedure may leave us with "kinks" in the itheta0 
	// variables. It may be that we need to recalculate those for all of the 
	// new points and the ones after them by making one more pass. I'm hoping
	// it is a realitively small correction though so we can skip it here.
	return 0;
}

//---------------------------------
// DistToRT
//---------------------------------
double DReferenceTrajectory::DistToRT(DVector3 hit)
{
if(Nswim_steps<1)_DBG__;
	// First, find closest step to point
	swim_step_t *swim_step = swim_steps;
	swim_step_t *step=NULL;
	double min_delta2 = 1.0E6;
	for(int i=0; i<Nswim_steps; i++, swim_step++){

		DVector3 pos_diff = swim_step->origin - hit;
		double delta2 = pos_diff.Mag2();
		if(delta2 < min_delta2){
			min_delta2 = delta2;
			step = swim_step;
		}
	}
	if(step==NULL){
		_DBG_<<"\"hit\" passed to DistToRT(DVector3) out of range!"<<endl;
		hit.Print();
		return 1.0E6;
	}
	
	// store last step
	last_swim_step=step;

	
	// Next, define a point on the helical segment defined by the
	// swim step it the RT coordinate system. The directions of
	// the RT coordinate system are defined by step->xdir, step->ydir,
	// and step->zdir. The coordinates of a point on the helix
	// in this coordinate system are:
	//
	//   x = Ro*(cos(phi) - 1)
	//   y = Ro*sin(phi)
	//   z = phi*(dz/dphi)
	//
	// where phi is the phi angle of the point in this coordinate system.
	// phi=0 corresponds to the swim step point itself
	//
	// Transform the given coordinates to the RT coordinate system
	// and call these x0,y0,z0. Then, the distance of point to a
	// point on the helical segment is given by:
	//
	//   d^2 = (x0-x)^2 + (y0-y)^2 + (z0-z)^2
	//
	// where x,y,z are all functions of phi as given above.
	//
	// writing out d^2 in terms of phi, but using the small angle
	// approximation for the trig functions, an equation for the
	// distance in only phi is obtained. Taking the derivative 
	// and setting it equal to zero leaves a 3rd order polynomial
	// in phi whose root corresponds to the minimum distance.
	// Skipping some math, this equation has the form:
	//
	// d(d^2)/dphi = 0 = Ro^2*phi^3 + 2*alpha*phi + beta
	//
	// where:
	//       alpha = x0*Ro + Ro^2 + (dz/dphi)^2
	//
	//        beta = -2*y0*Ro - 2*z0*(dz/dphi)
	//
	// The above 3rd order poly is convenient in that it does not
	// contain a phi^2 term. This means we can skip the step
	// done in the general case where a change of variables is
	// made such that the 2nd order term disappears.
	//
	// In general, an equation of the form
	//
	//  w^3 + 3.0*b*w + 2*c = 0 
	//
	// has one real root:
	//
	//  w0 = q - p
	//
	// where:
	//    q^3 = d - c
	//    p^3 = d + c
	//    d^2 = b^3 + c^2      (don't confuse with d^2 above!)
	//
	// So for us ...
	//
	//    3b = 2*alpha/(Ro^2)
	//    2c = beta/(Ro^2)

	hit -= step->origin;
	double x0 = hit.Dot(step->sdir);
	double y0 = hit.Dot(step->tdir);
	double z0 = hit.Dot(step->udir);
	double &Ro = step->Ro;
	double Ro2 = Ro*Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	double alpha = x0*Ro + Ro2 + pow(dz_dphi,2.0);
	double beta = -2.0*y0*Ro - 2.0*z0*dz_dphi;
	double b = (2.0*alpha/Ro2)/3.0;
	double c = (beta/Ro2)/2.0;
	double d = sqrt(pow(b,3.0) + pow(c,2.0));
	double q = pow(d-c, 1.0/3.0);
	double p = pow(d+c, 1.0/3.0);
	double phi = q - p;
	
	double dist2 = Ro2/4.0*pow(phi,4.0) + alpha*pow(phi,2.0) + beta*phi + x0*x0 + y0*y0 + z0*z0;

	return sqrt(dist2);
}

//---------------------------------
// FindClosestSwimStep
//---------------------------------
DReferenceTrajectory::swim_step_t* DReferenceTrajectory::FindClosestSwimStep(const DCoordinateSystem *wire, int *istep_ptr)
{
	/// Find the closest swim step to the given wire. The value of
	/// "L" should be the active wire length. The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	
	if(istep_ptr)*istep_ptr=-1;
	
	if(Nswim_steps<1){
		_DBG_<<"No swim steps! You must \"Swim\" the track before calling FindClosestSwimStep(...)"<<endl;
	}

	// Make sure we have a wire first!
	if(!wire)return NULL;
	
	// Loop over swim steps and find the one closest to the wire
	swim_step_t *swim_step = swim_steps;
	swim_step_t *step=NULL;
	double min_delta2 = 1.0E6;
	double L_over_2 = wire->L/2.0; // half-length of wire in cm
	int istep=-1;
	for(int i=0; i<Nswim_steps; i++, swim_step++){
		// Find the point's position along the wire. If the point
		// is past the end of the wire, calculate the distance
		// from the end of the wire.
		DVector3 pos_diff = swim_step->origin - wire->origin;
		double u = wire->udir.Dot(pos_diff);

		// Find distance perpendicular to wire
		double delta2 = pos_diff.Mag2() - u*u;

		// If point is past end of wire, calculate distance
		// from wire's end by adding on distance along wire direction.
		if(fabs(u)>L_over_2){
			delta2 += pow(fabs(u)-L_over_2, 2.0);
		}

		if(delta2 < min_delta2){
			min_delta2 = delta2;
			step = swim_step;
			istep=i;
		}
	}

	if(istep_ptr)*istep_ptr=istep;

	return step;	
}

//---------------------------------
// DistToRT
//---------------------------------
double DReferenceTrajectory::DistToRT(const DCoordinateSystem *wire, double *s)
{
	/// Find the closest distance to the given wire in cm. The value of
	/// "L" should be the active wire length (in cm). The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	swim_step_t *step=FindClosestSwimStep(wire);

	return step ? DistToRT(wire, step, s):std::numeric_limits<double>::quiet_NaN();
}

//---------------------------------
// DistToRTBruteForce
//---------------------------------
double DReferenceTrajectory::DistToRTBruteForce(const DCoordinateSystem *wire, double *s)
{
	/// Find the closest distance to the given wire in cm. The value of
	/// "L" should be the active wire length (in cm). The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	swim_step_t *step=FindClosestSwimStep(wire);
	
	return step ? DistToRTBruteForce(wire, step, s):std::numeric_limits<double>::quiet_NaN();
}

//------------------
// DistToRT
//------------------
double DReferenceTrajectory::DistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double *s)
{
	/// Calculate the distance of the given wire(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim step
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the wire.
	/// IMPORTANT: This approximates the helix locally by a parabola.
	/// This means the swim step should be fairly close
	/// to the wire so that this approximation is valid. If the
	/// reference trajectory from which the swim step came is too
	/// sparse, the results will not be nearly as good.
	
	// Interestingly enough, this is one of the harder things to figure
	// out in the tracking code which is why the explanations may be
	// a bit long.

	// The general idea is to define the helix in a coordinate system
	// in which the wire runs along the z-axis. The distance to the
	// wire is then defined just in the X/Y plane of this coord. system.
	// The distance is expressed as a function of the phi angle in the
	// natural coordinate system of the helix. This way, phi=0 corresponds
	// to the swim step point itself and the DOCA point should be
	// at a small phi angle.
	//
	// The minimum distance between the helical segment and the wire
	// will be a function of sin(phi), cos(phi) and phi. Approximating
	// sin(phi) by phi and cos(phi) by (1-phi^2) leaves a 4th order
	// polynomial in phi. Taking the derivative leaves a 3rd order
	// polynomial whose root is the phi corresponding to the 
	// Distance Of Closest Approach(DOCA) point on the helix. Plugging
	// that value of phi back into the distance formula gives
	// us the minimum distance between the track and the wire.

	// First, we need to define the coordinate system in which the 
	// wire runs along the z-axis. This is actually done already
	// in the CDC package for each wire once, at program start.
	// The directions of the axes are defined in wire->sdir,
	// wire->tdir, and wire->udir.
	
	// Next, define a point on the helical segment defined by the
	// swim step it the RT coordinate system. The directions of
	// the RT coordinate system are defined by step->xdir, step->ydir,
	// and step->zdir. The coordinates of a point on the helix
	// in this coordinate system are:
	//
	//   x = Ro*(cos(phi) - 1)
	//   y = Ro*sin(phi)
	//   z = phi*(dz/dphi)
	//
	// where phi is the phi angle of the point in this coordinate system.
	
	// Now, a vector describing the helical point in the LAB coordinate
	// system is:
	//
	//  h = x*xdir + y*ydir + z*zdir + pos
	//
	// where h,xdir,ydir,zdir and pos are all 3-vectors.
	// xdir,ydir,zdir are unit vectors defining the directions
	// of the RT coord. system axes in the lab coord. system.
	// pos is a vector defining the position of the swim step
	// in the lab coord.system 
	
	// Now we just need to find the extent of "h" in the wire's
	// coordinate system (period . means dot product):
	//
	// s = (h-wpos).sdir
	// t = (h-wpos).tdir
	// u = (h-wpos).udir
	//
	// where wpos is the position of the center of the wire in
	// the lab coord. system and is given by wire->wpos.

	// At this point, the values of s,t, and u repesent a point
	// on the helix in the coord. system of the wire with the
	// wire in the "u" direction and positioned at the origin.
	// The distance(squared) from the wire to the point on the helix
	// is given by:
	//
	// d^2 = s^2 + t^2
	//
	// where s and t are both functions of phi.

	// So, we'll define the values of "s" and "t" above as:
	//
	//   s = A*x + B*y + C*z + D
	//   t = E*x + F*y + G*z + H
	//
	// where A,B,C,D,E,F,G, and H are constants defined below
	// and x,y,z are all functions of phi defined above.
	// (period . means dot product)
	//
	// A = sdir.xdir
	// B = sdir.ydir
	// C = sdir.zdir
	// D = sdir.(pos-wpos)
	//
	// E = tdir.xdir
	// F = tdir.ydir
	// G = tdir.zdir
	// H = tdir.(pos-wpos)
	const DVector3 &xdir = step->sdir;
	const DVector3 &ydir = step->tdir;
	const DVector3 &zdir = step->udir;
	const DVector3 &sdir = wire->sdir;
	const DVector3 &tdir = wire->tdir;
	const DVector3 &udir = wire->udir;
	DVector3 pos_diff = step->origin - wire->origin;
	
	double A = sdir.Dot(xdir);
	double B = sdir.Dot(ydir);
	double C = sdir.Dot(zdir);
	double D = sdir.Dot(pos_diff);

	double E = tdir.Dot(xdir);
	double F = tdir.Dot(ydir);
	double G = tdir.Dot(zdir);
	double H = tdir.Dot(pos_diff);

	// OK, here is the dirty part. Using the approximations given above
	// to write the x and y functions in terms of phi^2 and phi (instead
	// of cos and sin) we put them into the equations for s and t above.
	// Then, inserting those into the equation for d^2 above that, we
	// get a very long equation in terms of the constants A,...H and
	// phi up to 4th order. Combining coefficients for similar powers
	// of phi yields an equation of the form:
	//
	// d^2 = Q*phi^4 + R*phi^3 + S*phi^2 + T*phi + U
	//
	// The dirty part is that it takes the better part of a sheet of
	// paper to work out the relations for Q,...U in terms of
	// A,...H, and Ro, dz/dphi. You can work it out yourself on
	// paper to verify that the equations below are correct.
	double Ro = step->Ro;
	double Ro2 = Ro*Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	double Q = pow(A*Ro/2.0, 2.0) + pow(E*Ro/2.0, 2.0);
	double R = -(2.0*A*B*Ro2 + 2.0*A*C*Ro*dz_dphi + 2.0*E*F*Ro2 + 2.0*E*G*Ro*dz_dphi)/2.0;
	double S = pow(B*Ro, 2.0) + pow(C*dz_dphi,2.0) + 2.0*B*C*Ro*dz_dphi - 2.0*A*D*Ro/2.0
					+ pow(F*Ro, 2.0) + pow(G*dz_dphi,2.0) + 2.0*F*G*Ro*dz_dphi - 2.0*E*H*Ro/2.0;
	double T = 2.0*B*D*Ro + 2.0*C*D*dz_dphi + 2.0*F*H*Ro + 2.0*G*H*dz_dphi;
	double U = D*D + H*H;
	
	// Aaarghh! my fingers hurt just from typing all of that!
	//
	// OK, now we differentiate the above equation for d^2 to get:
	//
	// d(d^2)/dphi = 4*Q*phi^3 + 3*R*phi^2 + 2*S*phi + T
	//
	// NOTE: don't confuse "R" with "Ro" in the above equations!
	//
	// Now we have to solve the 3rd order polynomial for the phi value of
	// the point of closest approach on the RT. This is a well documented
	// procedure. Essentially, when you have an equation of the form:
	//
	//  x^3 + a2*x^2 + a1*x + a0 = 0;
	//
	// a change of variables is made such that w = x + a2/3 which leads
	// to a third order poly with no w^2 term:
	//
	//  w^3 + 3.0*b*w + 2*c = 0 
	//
	// where:
	//    b = a1/3 - (a2^2)/9
	//    c = a0/2 - a1*a2/6  + (a2^3)/27
	//
	// The one real root of this is:
	//
	//  w0 = q - p
	//
	// where:
	//    q^3 = d - c
	//    p^3 = d + c
	//    d^2 = b^3 + c^2      (don't confuse with d^2 above!)
	//
	// For us this means that:
	//    a2 = 3*R/(4*Q)
	//    a1 = 2*S/(4*Q)
	//    a0 =   T/(4*Q)
	//
	// A potential problem could occur if Q is at or very close to zero.
	// This situation occurs when both A and E are zero. This would mean
	// that both sdir and tdir are perpendicular to xdir which means
	// xdir is in the same direction as udir (got that?). Physically,
	// this corresponds to the situation when both the momentum and
	// the magnetic field are perpendicular to the wire (though not
	// necessarily perpendicular to each other). This situation can't
	// really occur in the CDC detector where the chambers are well
	// contained in a region where the field is essentially along z as
	// are the wires.
	//
	// Just to be safe, we check that Q is greater than
	// some minimum before solving for phi. If it is too small, we fall
	// back to solving the quadratic equation for phi.
	double phi =0.0;
	if(fabs(Q)>1.0E-6){
		double fourQ = 4.0*Q;
		double a2 = 3.0*R/fourQ;
		double a1 = 2.0*S/fourQ;
		double a0 =     T/fourQ;

		double b = a1/3.0 - a2*a2/9.0;
		double c = a0/2.0 - a1*a2/6.0 + a2*a2*a2/27.0;

		double d = sqrt(pow(b, 3.0) + pow(c, 2.0)); // occasionally, this is zero. See below
		double q = pow(d - c, 1.0/3.0);
		double p = pow(d + c, 1.0/3.0);

		double w0 = q - p;
		phi = w0 - a2/3.0;
	}else{
		double a = 3.0*R;
		double b = 2.0*S;
		double c = 1.0*T;
		phi = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a); 
	}
	
	// The accuracy of this method is limited by how close the step is to the
	// actual minimum. If the value of phi is large then the step size is
	// not too close and we should add another couple of steps in the right
	// place in order to get a more accurate value. Note that while this will
	// increase the time it takes this round, presumably the fitter will be
	// calling this often for each wire and having a high density of points
	// near the wires will just make subsequent calls go quicker. This also
	// allows larger initial step sizes with the high density regions getting
	// filled in as needed leading to overall faster tracking.
	if(finite(phi) && fabs(phi)>2.0E-4){
		if(dist_to_rt_depth>=3){
			_DBG_<<"3 or more recursive calls to DistToRT(). Something is wrong! bailing ..."<<endl;
			//for(int k=0; k<Nswim_steps; k++){
			//	DVector3 &v = swim_steps[k].origin;
			//	_DBG_<<"  "<<k<<": "<<v.X()<<", "<<v.Y()<<", "<<v.Z()<<endl;
			//}
			//exit(-1);
			return std::numeric_limits<double>::quiet_NaN();
		}
		double scale_step = 1.0;
		double s_range = 1.0*scale_step;
		double step_size = 0.02*scale_step;
		int err = InsertSteps(step, phi>0.0 ? +s_range:-s_range, step_size);			// Add new steps near this step by swimming in the direction of phi
		if(!err){
			step=FindClosestSwimStep(wire);									// Find the new closest step
			if(!step)return std::numeric_limits<double>::quiet_NaN();
			dist_to_rt_depth++;
			double doca = DistToRT(wire, step, s);								// re-call ourself with the new step
			dist_to_rt_depth--;
			return doca;
		}else{
			if(err<0)return std::numeric_limits<double>::quiet_NaN();
			
			// If InsertSteps() returns an error > 0 then it indicates that it
			// was unable to add additional steps (perhaps because there
			// aren't enough spaces available). In that case, we just go ahead
			// and use the phi we have and make the best estimate possible.
		}
	}
	
	// It is possible at this point that the value of phi corresponds to
	// a point past the end of the wire. We should check for this here and
	// recalculate, if necessary, the DOCA at the end of the wire. First,
	// calculate h (the vector defined way up above) and dot it into the
	// wire's u-direction to get the position of the DOCA point along the
	// wire.
	double x = -Ro*phi*phi/2.0;
	double y = Ro*phi;
	double z = dz_dphi*phi;
	DVector3 h = pos_diff + x*xdir + y*ydir + z*zdir;
	double u = h.Dot(udir);
	if(fabs(u) > wire->L/2.0){
		// Looks like our DOCA point is past the end of the wire.
		// Find phi corresponding to the end of the wire.
		double L_over_2 = u>0.0 ? wire->L/2.0:-wire->L/2.0;
		double a = -Ro*udir.Dot(xdir)/2.0;
		double b = Ro*udir.Dot(ydir) + dz_dphi*udir.Dot(zdir);
		double c = udir.Dot(pos_diff) - L_over_2;
		double phi1 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a); 
		double phi2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
		phi = fabs(phi1)<fabs(phi2) ? phi1:phi2;
		u=L_over_2;
	}
	this->last_dist_along_wire = u;

	// Sometimes the "d" used in solving the 3rd order polynmial above 
	// can be nan due to the sqrt argument being negative. I'm not sure
	// exactly what this means (e.g. is this due to round-off, are there
	// no roots, ...) Also unclear is how to handle it. The only two choices
	// I can think of are : 1.) set phi to zero or 2.) return the nan
	// value. Option 1.) tries to keep the hit while option 2 ties to ignore
	// it. Both options should probably be studied at some point. For now
	// though, it looks (at least preliminarily) like this occurs slightly
	// less than 1% of the time on valid hits so we go ahead with option 2.

	// Use phi to calculate DOCA
	double d2 = U + phi*(T + phi*(S + phi*(R + phi*Q)));
	double d = sqrt(d2);
	
	// Calculate distance along track ("s")
	double dz = dz_dphi*phi;
	double Rodphi = Ro*phi;
	double ds = sqrt(dz*dz + Rodphi*Rodphi);
	if(s)*s=step->s + (phi>0.0 ? ds:-ds);
	
	// Remember phi and step so additional info on the point can be obtained
	this->last_phi = phi;
	this->last_swim_step = step;
	this->last_dz_dphi = dz_dphi;

	return d; // WARNING: This could return nan!
}

//------------------
// DistToRTBruteForce
//------------------
double DReferenceTrajectory::DistToRTBruteForce(const DCoordinateSystem *wire, const swim_step_t *step, double *s)
{
	/// Calculate the distance of the given wire(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim step
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the wire.
	/// IMPORTANT: This calculates the distance using a "brute force"
	/// method of taking tiny swim steps to find the minimum distance.
	/// It is vey SLOW and you should be using DistToRT(...) instead.
	/// This is only here to provide an independent check of DistToRT(...).
	
	const DVector3 &xdir = step->sdir;
	const DVector3 &ydir = step->tdir;
	const DVector3 &zdir = step->udir;
	const DVector3 &sdir = wire->sdir;
	const DVector3 &tdir = wire->tdir;
	DVector3 pos_diff = step->origin - wire->origin;
	
	double Ro = step->Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	// Brute force
	double min_d2 = 1.0E6;
	double phi=M_PI;
	for(int i=-2000; i<2000; i++){
		double myphi=(double)i*0.000005;
		DVector3 d = Ro*(cos(myphi)-1.0)*xdir
	            	+ Ro*sin(myphi)*ydir
						+ dz_dphi*myphi*zdir
						+ pos_diff;

		double d2 = pow(d.Dot(sdir),2.0) + pow(d.Dot(tdir),2.0);
		if(d2<min_d2){
			min_d2 = d2;
			phi = myphi;
			this->last_phi = myphi;
		}
	}
	double d2 = min_d2;
	double d = sqrt(d2);
	this->last_swim_step = step;
	this->last_dz_dphi = dz_dphi;
	
	// Calculate distance along track ("s")
	double dz = dz_dphi*phi;
	double Rodphi = Ro*phi;
	double ds = sqrt(dz*dz + Rodphi*Rodphi);
	if(s)*s=step->s + (phi>0.0 ? ds:-ds);

	return d;
}

//------------------
// Straw_dx
//------------------
double DReferenceTrajectory::Straw_dx(const DCoordinateSystem *wire, double radius)
{
	/// Find the distance traveled within the specified radius of the
	/// specified wire. This will give the "dx" component of a dE/dx
	/// measurement for cylindrical geometry as we have with straw tubes.
	///
	/// At this point, the estimate is done using a simple linear 
	/// extrapolation from the DOCA point in the direction of the momentum
	/// to the 2 points at which it itersects the given radius. Segments
	/// which extend past the end of the wire will be clipped to the end
	/// of the wire before calculating the total dx.
	
	// First, find the DOCA point for this wire
	double s;
	double doca = DistToRT(wire, &s);
	
	// If doca is outside of the given radius, then we're done
	if(doca>=radius)return 0.0;
	
	// Get the location and momentum direction of the DOCA point
	DVector3 pos, momdir;
	GetLastDOCAPoint(pos, momdir);
	momdir.SetMag(1.0);
	
	// Get wire direction
	const DVector3 &udir = wire->udir;
	
	// Calculate vectors used to form quadratic equation for "alpha"
	// the distance along the mometum direction from the DOCA point
	// to the intersection with a cylinder of the given radius.
	DVector3 A = udir.Cross(pos-wire->origin);
	DVector3 B = udir.Cross(momdir);
	
	// If the magnitude of B is zero at this point, it means the momentum
	// direction is parallel to the wire. In this case, this method will
	// not work. Return NaN.
	if(B.Mag()<1.0E-10)return std::numeric_limits<double>::quiet_NaN();
	
	double a = B.Mag();
	double b = A.Dot(B);
	double c = A.Mag() - radius;
	double d = sqrt(b*b - 4.0*a*c);
	
	// The 2 roots should correspond to the 2 intersection points.
	double alpha1 = (-b + d)/(2.0*a);
	double alpha2 = (-b - d)/(2.0*a);
	
	DVector3 int1 = pos + alpha1*momdir;
	DVector3 int2 = pos + alpha2*momdir;
	
	// Check if point1 is past the end of the wire
	double q = udir.Dot(int1 - wire->origin);
	if(fabs(q) > wire->L/2.0){
		double gamma = udir.Dot(wire->origin - pos) + (q>0.0 ? +1.0:-1.0)*wire->L/2.0;
		gamma /= momdir.Dot(udir);
		int1 = pos + gamma*momdir;
	}

	// Check if point2 is past the end of the wire
	q = udir.Dot(int2 - wire->origin);
	if(fabs(q) > wire->L/2.0){
		double gamma = udir.Dot(wire->origin - pos) + (q>0.0 ? +1.0:-1.0)*wire->L/2.0;
		gamma /= momdir.Dot(udir);
		int2 = pos + gamma*momdir;
	}
	
	// Calculate distance
	DVector3 delta = int1 - int2;
	
	return delta.Mag();
}

//------------------
// GetLastDOCAPoint
//------------------
void DReferenceTrajectory::GetLastDOCAPoint(DVector3 &pos, DVector3 &mom)
{
	/// Use values saved by the last call to one of the DistToRT functions
	/// to calculate the 3-D DOCA position in lab coordinates and momentum
	/// in GeV/c.
	const DVector3 &xdir = last_swim_step->sdir;
	const DVector3 &ydir = last_swim_step->tdir;
	const DVector3 &zdir = last_swim_step->udir;

	double x = -(last_swim_step->Ro/2.0)*last_phi*last_phi;
	double y = last_swim_step->Ro*last_phi;
	double z = last_dz_dphi*last_phi;

	pos = last_swim_step->origin + x*xdir + y*ydir + z*zdir;
	mom = last_swim_step->mom;
	mom.Rotate(-last_phi, zdir);
}

//------------------
// GetLastDOCAPoint
//------------------
DVector3 DReferenceTrajectory::GetLastDOCAPoint(void)
{
	/// Use values saved by the last call to one of the DistToRT functions
	/// to calculate the 3-D DOCA position in lab coordinates. This is
	/// mainly intended for debugging.
	const DVector3 &xdir = last_swim_step->sdir;
	const DVector3 &ydir = last_swim_step->tdir;
	const DVector3 &zdir = last_swim_step->udir;
	double Ro = last_swim_step->Ro;
	double delta_z = last_swim_step->mom.Dot(zdir);
	double delta_phi = last_swim_step->mom.Dot(ydir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	double x = -(Ro/2.0)*last_phi*last_phi;
	double y = Ro*last_phi;
	double z = dz_dphi*last_phi;

	return last_swim_step->origin + x*xdir + y*ydir + z*zdir;
}
