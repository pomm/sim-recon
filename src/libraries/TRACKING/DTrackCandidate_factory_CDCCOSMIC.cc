// $Id$
//
//    File: DTrackCandidate_factory_CDCCOSMIC.cc
// Created: Sat Jun 28 16:50:07 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.2.0 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include "Math/Minimizer.h"
#include "TMinuitMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
using namespace std;

#include <CDC/DCDCTrackHit.h>
#include <TRACKING/DTrackFinder.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALShower.h>
#include "DTrackCandidate_factory_CDCCOSMIC.h"
using namespace jana;

// This simple track finder is made for fitting straight line cosmic tracks 
// in the CDC. It will automatically make a track from all CDC hits it finds. 
// It does this by doing a linear regression on the CDC axial hits
// in the X/Y plane.
// Updated to include timing information in the fit.
// Will use TMinuit (Not recommended for track fitting, but this isn't really tne "production" method)

bool DTrackCandidate_CDCCOSMIC_cdc_hit_cmp(const DCDCTrackHit *a,
        const DCDCTrackHit *b){

    return(a->wire->origin.Y()>b->wire->origin.Y());
}

vector<double> Measurements, MeasurementErrors;
vector<const DCDCWire *> Wires;

Double_t fit_function(const DCDCWire * wire,const Double_t *par)
{
    // par[4]={x0, y0, xslope, yslope}
    double udirx = wire->udir.x();
    double udiry = wire->udir.y();
    double udirz = wire->udir.z();
    double wx = wire->origin.x();
    double wy = wire->origin.y();
    double wz = wire->origin.z();
    // Distance of closest approach from track to wire
    double value = sqrt(pow(udiry*(wx - par[0] - wz*par[2]) + udirx*(-wy + par[1] + wz*par[3]) + udirz*(wy*par[2] - par[2]*par[1] - wx*par[3] + par[0]*par[3]),2)/
            (pow(udiry,2)*(1 + pow(par[2],2)) - 2*udiry*udirz*par[3] - 2*udirx*par[2]*(udirz + udiry*par[3]) + 
             pow(udirx,2)*(1 + pow(par[3],2)) + pow(udirz,2)*(pow(par[2],2) + pow(par[3],2))));
    return value;
}

double calc_chi_square(const Double_t *par)
{
    //calculate chisquare
    double chisq = 0;
    for (unsigned int i=0;i<Measurements.size(); i++) {
        // chi square is the quadratic sum of the distance from the point to the function weighted by its error
        double delta  = (Measurements[i]-fit_function(Wires[i],par))/MeasurementErrors[i];
        chisq += delta*delta;
    }
    return chisq;
}

// Convert time to distance for the cdc
double DTrackCandidate_factory_CDCCOSMIC::CDCDriftDistance(double t){
    double d=0.;
    /*
    // Try 9th order polynomial fit to the drift distribution
    double p[10]= { 0.0, 0.0077461, -9.31638e-05, 7.3452e-07, -3.47947e-09, 1.01626e-11, -1.84211e-14, 2.01549e-17, -1.21754e-20, 3.11542e-24};
    if (t > 0.0){
        d = p[0] + p[1] * t + p[2] * pow(t,2) + p[3] * pow(t,3) + p[4] * pow(t,4)
            + p[5] * pow(t,5) + p[6] * pow(t,6) + p[7] * pow(t,7) + p[8] * pow(t,8) + p[9] * pow(t,9);
    }
    */
    
    if (t>cdc_drift_table[cdc_drift_table.size()-1]) return 0.78;
    if (t>0){
        unsigned int index=0;
        index=Locate(cdc_drift_table,t);
        double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
        double frac=(t-cdc_drift_table[index])/dt;
        d=0.01*(double(index)+frac); 
    }
    return d;
}

// Smearing function derived from fitting residuals
inline double DTrackCandidate_factory_CDCCOSMIC::CDCDriftVariance(double t){ 
    //double sigma = 0.015;
    //return sigma*sigma;
    //  return 0.001*0.001;
    //if (t<0.) t=0.;
    //cout << "Entering Variance lookup for time " << t << endl;

    double cutoffTime = 5.0;
    double V = 0.0507;
    if (t>0){
        if (t>cdc_drift_table_max){
            // Here the drift time is too large. Return the edge of the straw with a large error
            V=0.0507; // straw radius^2 / 12
            return V;
        }
        //cout << "The resolution parameters are {" << CDC_RES_PAR1 <<","<< CDC_RES_PAR2 << "}" << endl;
        double sigma=CDC_RES_PAR1/(t+1.)+CDC_RES_PAR2;
        //cout << "Sigma = " << sigma << endl;
        // This function is very poorly behaved at low drift times
        // For times less than cutoffTime assume a linear behavior of our function.
       // if( t < cutoffTime ){
       //     double slope = -1.0 * CDC_RES_PAR1 / (( cutoffTime + 1) * (cutoffTime + 1));
       //     sigma = (CDC_RES_PAR1/(cutoffTime+1.)+CDC_RES_PAR2) + slope * (t - cutoffTime);
            //cout << "Earlier than cutoff time...Sigma = " << sigma << endl;
        //}

        V=sigma*sigma;
        return V;
    }
    else { // Time is negative, or exactly zero, choose position at wire, with error of t=0 hit
        double slope = -1.0 * CDC_RES_PAR1 / (( cutoffTime + 1) * (cutoffTime + 1));
        double sigma = (CDC_RES_PAR1/(cutoffTime+1.)+CDC_RES_PAR2) + slope * (0.0 - cutoffTime);
        //cout << "Time is negative...sigma = " << sigma << endl;
        V=sigma*sigma; //Should include T0 variance...
        return V;
        //V=0.0507; // straw radius^2 / 12
    }
    //cout << "Somehow we got here...returning Variance = " << V << endl;
    return V;
    
}

double DTrackCandidate_factory_CDCCOSMIC::CDCTrackError(const DCDCWire *wire , const double *xs, double *dxs){
    //First the wire parameters
    double udirx = wire->udir.x();
    double udiry = wire->udir.y();
    double udirz = wire->udir.z();
    double wx = wire->origin.x();
    double wy = wire->origin.y();
    double wz = wire->origin.z();

    //The fitted track parameters
    double x0 = xs[0];
    double y0 = xs[1];
    double xslope = xs[2];
    double yslope = xs[3];

    // We need the jacobian "matrix" for this measurement
    // Taking these directly form Mathematica
    
    //df/dx0
    double dfdx0 = ((udiry - udirz*yslope)*sqrt(pow(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + 
                             udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope),2)/
                       (pow(udiry,2)*(1 + pow(xslope,2)) - 2*udiry*udirz*yslope - 2*udirx*xslope*(udirz + udiry*yslope) + pow(udirx,2)*(1 + pow(yslope,2)) + 
                                 pow(udirz,2)*(pow(xslope,2) + pow(yslope,2)))))/
           (udiry*(-wx + x0 + wz*xslope) + (udirx - udirz*xslope)*(wy - y0) - (udirx*wz + udirz*(-wx + x0))*yslope);

    //df/dy0
    double dfdy0 = ((udirx - udirz*xslope)*sqrt(pow(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + 
                             udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope),2)/
                       (pow(udiry,2)*(1 + pow(xslope,2)) - 2*udiry*udirz*yslope - 2*udirx*xslope*(udirz + udiry*yslope) + pow(udirx,2)*(1 + pow(yslope,2)) + 
                                 pow(udirz,2)*(pow(xslope,2) + pow(yslope,2)))))/
           (udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope));

    //df/dxslope
    double dfdxslope = ((udiry - udirz*yslope)*pow(pow(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + 
                             udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope),2)/
                       (pow(udiry,2)*(1 + pow(xslope,2)) - 2*udiry*udirz*yslope - 2*udirx*xslope*(udirz + udiry*yslope) + pow(udirx,2)*(1 + pow(yslope,2)) + 
                                 pow(udirz,2)*(pow(xslope,2) + pow(yslope,2))),1.5)*
                 (pow(udiry,2)*(wz + wx*xslope - x0*xslope) - udiry*udirz*(wy - y0 + wz*yslope) + pow(udirx,2)*(wz + (wy - y0)*yslope) + 
                         pow(udirz,2)*(wx*xslope - x0*xslope + wy*yslope - y0*yslope) - 
                                udirx*(udirz*(wx - x0 + wz*xslope) + udiry*(wy*xslope - xslope*y0 + wx*yslope - x0*yslope))))/
           ((udiry*(-wx + x0 + wz*xslope) + (udirx - udirz*xslope)*(wy - y0) - (udirx*wz + udirz*(-wx + x0))*yslope)*
                 pow(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope),2));

    //df/dyslope
    double dfdyslope = ((udirx - udirz*xslope)*(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope))*
                 (pow(udiry,2)*(wz + wx*xslope - x0*xslope) - udiry*udirz*(wy - y0 + wz*yslope) + pow(udirx,2)*(wz + (wy - y0)*yslope) + 
                         pow(udirz,2)*(wx*xslope - x0*xslope + wy*yslope - y0*yslope) - 
                                udirx*(udirz*(wx - x0 + wz*xslope) + udiry*(wy*xslope - xslope*y0 + wx*yslope - x0*yslope))))/
           (sqrt(pow(udiry*(wx - x0 - wz*xslope) + udirx*(-wy + y0 + wz*yslope) + udirz*(wy*xslope - xslope*y0 - wx*yslope + x0*yslope),2)/
                        (pow(udiry,2)*(1 + pow(xslope,2)) - 2*udiry*udirz*yslope - 2*udirx*xslope*(udirz + udiry*yslope) + pow(udirx,2)*(1 + pow(yslope,2)) + 
                                  pow(udirz,2)*(pow(xslope,2) + pow(yslope,2))))*
                 pow(pow(udiry,2)*(1 + pow(xslope,2)) - 2*udiry*udirz*yslope - 2*udirx*xslope*(udirz + udiry*yslope) + 
                            pow(udirx,2)*(1 + pow(yslope,2)) + pow(udirz,2)*(pow(xslope,2) + pow(yslope,2)),2));

    // Since these only will be useful when we have the covariance matrix of the track parameters we have an aditional cut on our track quality.
    // Propogation of errors says that the error in the track measurement is given by
    // ( dfdx0 dfdy0 dfdxslope dfdyslope) (                         ) (  dfdx0  ) 
    //                                    ( Trk. covariance matrix  ) (  dfdy0  )  = The variance in the estimated track DOCA for this particular wire
    //                                    (                         ) (dfdxslope)
    //                                    (                         ) (dfdyslope)
    // If these were passed around as vectors and matrices this would be prettier, but since its short we can just type it out
    
    double variance = (dfdx0*dxs[0]  + dfdy0*dxs[1]  + dfdxslope*dxs[2]  + dfdyslope*dxs[3])  * dfdx0 
                    + (dfdx0*dxs[4]  + dfdy0*dxs[5]  + dfdxslope*dxs[6]  + dfdyslope*dxs[7])  * dfdy0
                    + (dfdx0*dxs[8]  + dfdy0*dxs[9]  + dfdxslope*dxs[10] + dfdyslope*dxs[11]) * dfdxslope
                    + (dfdx0*dxs[12] + dfdy0*dxs[13] + dfdxslope*dxs[14] + dfdyslope*dxs[15]) * dfdyslope;

    double error = sqrt(variance);

    return error;
}

void DTrackCandidate_factory_CDCCOSMIC::GetDOCAPhiandZ(const DCDCWire *wire, DTrackCandidate *locTrack, double &phi, double &z){
    // Get the vector pointing from the wire to the doca point
    DVector3 trackPosition = locTrack->position();
    DVector3 trackMomentum = locTrack->momentum();
    DVector3 wirePosition = wire->origin;
    DVector3 wireDirection = wire->udir;
    //wirePosition.Print(); wireDirection.Print();
    Float_t a = trackMomentum.Dot(trackMomentum);
    Float_t b = trackMomentum.Dot(wireDirection);
    Float_t c = wireDirection.Dot(wireDirection);
    DVector3 w0 = trackPosition - wirePosition;
    Float_t d = trackMomentum.Dot(w0);
    Float_t e = wireDirection.Dot(w0);
    Float_t sc = ((b*e - c*d)/(a*c-b*b));
    //if (sc < 0) continue; // Track must come from location away from origin
    DVector3 POCAOnTrack = trackPosition + sc * trackMomentum;
    z = POCAOnTrack.Z();
    DVector3 LOCA = w0 + ((b*e - c*d)/(a*c-b*b))*trackMomentum - ((a*e - b*d)/(a*c-b*b))*wireDirection;
    phi = LOCA.Phi();
    return;
}
// Locate a position in vector xx given x
unsigned int DTrackCandidate_factory_CDCCOSMIC::Locate(vector<double>&xx,
        double x){
    int n=xx.size();
    if (x==xx[0]) return 0;
    else if (x==xx[n-1]) return n-2;

    int jl=-1;
    int ju=n;
    int ascnd=(xx[n-1]>=xx[0]);
    while(ju-jl>1){
        int jm=(ju+jl)>>1;
        if ( (x>=xx[jm])==ascnd)
            jl=jm;
        else
            ju=jm;
    } 
    return jl;
}

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::init(void)
{
    bfield = new DMagneticFieldMapNoField(japp);
    rt = new DReferenceTrajectory(bfield);

    EXCLUDERING=0;
    if (gPARMS){
        gPARMS->SetDefaultParameter("CDCCOSMIC:EXCLUDERING", EXCLUDERING, "Ring to exclude in CDC Cosmic Tracking");
    }
    //rt->Rmax_interior = 100.0; // (cm)  set larger swim volume so we can swim through the BCAL 
    //rt->Rmax_exterior = 200.0; // (cm)  set larger swim volume so we can swim through the BCAL 

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::brun(jana::JEventLoop *eventLoop, int runnumber)
{
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    JCalibration *jcalib = dapp->GetJCalibration(runnumber);
    // Get pointer to TrackFinder object 
    vector<const DTrackFinder *> finders;
    eventLoop->Get(finders);

    if(finders.size()<1){
        _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
        return RESOURCE_UNAVAILABLE;
    }

    // Drop the const qualifier from the DTrackFinder pointer
    finder = const_cast<DTrackFinder*>(finders[0]);

    typedef map<string,double>::iterator iter_double;
    vector< map<string, double> > tvals;
    if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
        for(unsigned int i=0; i<tvals.size(); i++){
            map<string, double> &row = tvals[i];
            iter_double iter = row.find("t");
            cdc_drift_table.push_back(1000.*iter->second);
        }
    }

    if(cdc_drift_table.empty()){
        jerr << endl;
        jerr << " No values found for \"CDC/cdc_drift_table\"!" <<endl;
        jerr << endl;
        jerr << " This probably means you'r using an old calibration DB." << endl;
        jerr << " Check your JANA_CALIB_URL environment variable." << endl;
        jerr << " (This message printed from DCDCTrackHit_factory::brun())" << endl;
        exit(-1);
    }
    cdc_drift_table_min = cdc_drift_table[0];
    cdc_drift_table_max = cdc_drift_table[cdc_drift_table.size()-1];

    map<string, double> cdc_res_parms;
    jcalib->Get("CDC/cdc_resolution_parms", cdc_res_parms);
    CDC_RES_PAR1 = cdc_res_parms["res_par1"];
    CDC_RES_PAR2 = cdc_res_parms["res_par2"];

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::evnt(JEventLoop *loop, int eventnumber)
{

    vector <const DBCALShower *> bcalShowerVector;
    loop->Get(bcalShowerVector);

    bool hasBCALGuess = false;
    DMatrix4x1 BCALGuess;
    double maxEnergy = 0.;
    for (unsigned int i=0; i < bcalShowerVector.size(); i++){
        const DBCALShower *firstBCALShower = bcalShowerVector[i];
        // We need to form pairs of BCAL Showers to seed the tracks
        for (unsigned int j=i+1; j < bcalShowerVector.size(); j++){ 
            const DBCALShower *secondBCALShower = bcalShowerVector[j];
            // Some quality cuts on the BCAL Showers
            if (firstBCALShower->N_cell < 2 || secondBCALShower->N_cell < 2) continue;
            double x1 = firstBCALShower->x;
            double y1 = firstBCALShower->y;
            double z1 = firstBCALShower->z;
            double x2 = secondBCALShower->x;
            double y2 = secondBCALShower->y;
            double z2 = secondBCALShower->z;
            double xslope = (x2-x1) / (z2-z1);
            double yslope = (y2-y1) / (z2-z1);
            double x0 = x1 - xslope*z1;
            double y0 = y1 - yslope*z1;
            //cout << "BCAL Track parameter guess" << endl;
            //cout << "{x0,y0,xslope,yslope} = {" << x0 << "," << y0 << "," << xslope << "," << yslope << "}"<<endl;  
            double totalEnergy = firstBCALShower->E + secondBCALShower->E;
            if (totalEnergy > maxEnergy){
                DMatrix4x1 thisGuess(x0, y0, xslope,yslope);
                BCALGuess = thisGuess;
                hasBCALGuess=true;
            }
        }
    }

    // Use the track finder to get the tracks to attempt a fit
    // Look for tracks in the CDC
    vector<const DCDCTrackHit*>cdcs;
    loop->Get(cdcs);

    // Reset the track finder
    finder->Reset();

    if (cdcs.size()>4){
        for (size_t i=0;i<cdcs.size();i++) finder->AddHit(cdcs[i]);
        finder->FindAxialSegments();
        finder->LinkCDCSegments();

        const vector<DTrackFinder::cdc_track_t>tracks=finder->GetCDCTracks();
        if (tracks.size() != 1) return NOERROR;
        for (size_t i=0;i<tracks.size();i++){
            // Initial guess for state vector
            DMatrix4x1 S(tracks[i].S);
            // list of axial and stereo hits for this track
            hits=tracks[i].axial_hits;
            hits.insert(hits.end(),tracks[i].stereo_hits.begin(),
                    tracks[i].stereo_hits.end());
            sort(hits.begin(),hits.end(),DTrackCandidate_CDCCOSMIC_cdc_hit_cmp);

            // Use earliest cdc time to estimate t0
            double t0=1e6;
            for (unsigned int j=0;j<hits.size();j++){
                if (hits[j]->wire->ring == EXCLUDERING) continue;    
                double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
                double t_test=hits[j]->tdrift-L/29.98;
                if (t_test<t0) t0=t_test;
            }
            //Loop over hits
            Measurements.clear();
            MeasurementErrors.clear();
            Wires.clear();
            for (unsigned int j=0;j<hits.size();j++){
                //Calulate the Measurement from the drift time
                double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
                double tcorr = hits[j]->tdrift - L/29.98 - t0;
                double measurement = CDCDriftDistance(tcorr);
                Measurements.push_back(measurement);
                MeasurementErrors.push_back(sqrt(CDCDriftVariance(tcorr)));
                Wires.push_back(hits[j]->wire);
            }

            //Perform the fit 
            // create minimizer giving a name and a name (optionally) for the specific
            // algorithm
            // possible choices are: 
            //     minName                  algoName
            // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
            //  Minuit2                     Fumili2
            //  Fumili
            //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
            //                              BFGS2, SteepestDescent
            //  GSLMultiFit
            //   GSLSimAn
            //   Genetic
            ROOT::Math::Minimizer* min = 
                ROOT::Math::Factory::CreateMinimizer();
            if (min == NULL || min == 0){
                cout << "Error creating minimizer " << endl;
                return NOERROR;
            }
            //TMinuitMinimizer *min = new TMinuitMinimizer(); 

            // set tolerance , etc...
            min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
            //min->SetMaxIterations(10000);  // for GSL 
            min->SetTolerance(0.001);
            min->SetPrintLevel(-1);

            // create funciton wrapper for minmizer
            // a IMultiGenFunction type 
            ROOT::Math::Functor f(&calc_chi_square,4); 
            double step[4] = {0.01,0.01,0.01,0.01};
            // starting point

            double variable[4] = { S(0),S(1), S(2), S(3)};

            min->SetFunction(f);

            // Set the free variables to be minimized!
            min->SetVariable(0,"x0",variable[0], step[0]);
            min->SetVariable(1,"y0",variable[1], step[1]);
            min->SetVariable(2,"xslope",variable[2], step[2]);
            min->SetVariable(3,"yslope",variable[3], step[3]);

            // do the minimization
            if(!min->Minimize()) {
                // The fit failed if we got here
                // Check if there is a BCAL seed, and retry the fit 
                if (hasBCALGuess){
                    min->SetVariable(0,"x0",BCALGuess(0), step[0]);
                    min->SetVariable(1,"y0",BCALGuess(1), step[1]);
                    min->SetVariable(2,"xslope",BCALGuess(2), step[2]);
                    min->SetVariable(3,"yslope",BCALGuess(3), step[3]);
                    if(!min->Minimize()) return NOERROR; //Minimization Failed with BCAL Guess
                }    
                else return NOERROR;//Minimization Failed and no BCAL guess
            }

            // After the first pass at the fit, do some outlier rejection to improve the results
            float chi2cut = 20.0;

            const double *xs = min->X();

            //Loop over hits
            Measurements.clear();
            MeasurementErrors.clear();
            Wires.clear();
            vector<const DCDCTrackHit *> AcceptedHits;
            unsigned int nBeforeHits = hits.size(); 
            unsigned int nAfterHits = 0;
            for (unsigned int j=0;j<hits.size();j++){
                //Calulate the Measurement from the drift time
                double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
                double tcorr = hits[j]->tdrift - L/29.98 - t0;
                double measurement = CDCDriftDistance(tcorr);
                double residual = measurement - fit_function(hits[j]->wire,xs);
                double chi2calc = residual * residual / CDCDriftVariance(tcorr);
                if (chi2calc > chi2cut) continue;
                AcceptedHits.push_back(hits[j]);
                // Skip this hit if it needs to be exluded in the final fit
                // We still want it in the puls so mighat as well leave it in the accepted hits.
                if (hits[j]->wire->ring == EXCLUDERING) continue;
                nAfterHits++;
                Measurements.push_back(measurement);
                MeasurementErrors.push_back(sqrt(CDCDriftVariance(tcorr)));
                Wires.push_back(hits[j]->wire);
            }

            double acceptedFraction = (double) nAfterHits / nBeforeHits;
            if (acceptedFraction < 0.5 || nAfterHits < 5) return NOERROR;
            hits = AcceptedHits;

            // Refit
            // Set the free variables to be minimized!
            min->SetVariable(0,"x0",xs[0], step[0]);
            min->SetVariable(1,"y0",xs[1], step[1]);
            min->SetVariable(2,"xslope",xs[2], step[2]);
            min->SetVariable(3,"yslope",xs[3], step[3]);

            // do the minimization
            if(!min->Minimize() || min->CovMatrixStatus() != 3) { // Minimization converged and covariance matric accurate to get past here
                return NOERROR; //Minimization Failed after outlier rejection
            }

            xs = min->X();

            // Retrieve the covariance matrix
            double dxs[16];// covariance matrix stored as an array
            min->GetCovMatrix(dxs);

            DVector3 pos(xs[0], xs[1], 0);
            DVector3 mom(xs[2], xs[3], 1);

            DTrackCandidate *can = new DTrackCandidate();
            can->setMomentum(mom);
            can->setPosition(pos);
            can->setCharge(1.0); // make track draw as stright line
            can->setMass(0.139);
            can->setPID(PiPlus);
            can->chisq=min->MinValue();
            can->Ndof=nAfterHits - 4;
            can->rt = rt;
            rt->Swim(pos, mom, 1.0);

            // Loop through the hits and fill the pulls
            for (unsigned int j=0;j<hits.size();j++){
                // The following data is stored in the pulls
                // pull_t(double resi, double err,double s=0.0,
                //           double tdrift=0.0, double d=0.0,
                //           const DCDCTrackHit *cdc_hit=NULL,
                //           const DFDCPseudo *fdc_hit=NULL));
                double L=(hits[0]->wire->origin-hits[j]->wire->origin).Perp();
                double time = hits[j]->tdrift - L/29.98 - t0;
                double DOCA = fit_function(hits[j]->wire,xs);
                double measurement = CDCDriftDistance(time);
                double residual = measurement - DOCA;   
                double measurementError = sqrt(CDCDriftVariance(time));
                double trackError = CDCTrackError(hits[j]->wire, xs, dxs);
                double error = sqrt(pow(measurementError,2) + pow(trackError,2));
                //cout << " The error is as follows, measurment error = " << measurementError << " Tracking Error " << trackError << " Total " << error << endl;
                double docaphi, docaz;
                GetDOCAPhiandZ(hits[j]->wire, can, docaphi, docaz);
                can->pulls.push_back(DTrackFitter::pull_t(residual, error, 0.0, time , measurement, hits[j], NULL,docaphi,docaz));

            }
            _data.push_back(can);
        }
    }

    return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_CDCCOSMIC::fini(void)
{
    return NOERROR;
}

