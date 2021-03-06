// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <stdio.h>
#include <unistd.h>

#include "MyProcessor.h"
#include "PID/DPhoton.h"
#include "PID/DKinFit.h"
#include "TRACKING/DMCThrown.h"


int PAUSE_BETWEEN_EVENTS = 1;
int SKIP_BORING_EVENTS = 0;
int PRINT_ALL=0;
float MAX_EVENTS = 1e20;
int COUNT = 0;
//TFile *fout = NULL;
//TH1F *hpi0[4];

vector<string> toprint;

#define ansi_escape		((char)0x1b)
#define ansi_bold 		ansi_escape<<"[1m"
#define ansi_black		ansi_escape<<"[30m"
#define ansi_red			ansi_escape<<"[31m"
#define ansi_green		ansi_escape<<"[32m"
#define ansi_blue			ansi_escape<<"[34m"
#define ansi_normal		ansi_escape<<"[0m"
#define ansi_up(A)		ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)		ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)	ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)		ansi_escape<<"["<<(A)<<"D"

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
  char name[256];
  ///< Called once at program start.
  fout = new TFile(OUTNAME, "RECREATE","Output file");
  cout << "Opened ROOT file " << OUTNAME <<endl;
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<4;j++)
    {
      for(int k=0;k<5;k++)
      {
        sprintf(name,"hmm%d_%d_%d",i,j,k);
        hmm[i][j][k] = new TH1F(name,name,100,-0.3,0.3);

        sprintf(name,"hcl%d_%d_%d",i,j,k);
        hcl[i][j][k] = new TH1F(name,name,100,0.0,1.0);

        for(int p=0;p<4;p++)
        {
          sprintf(name,"hinvmass%d_%d_%d_%d",i,j,k,p);
          hinvmass[i][j][k][p] = new TH1F(name,name,100,0.3,0.7);
        }

        for(int p=0;p<24;p++)
        {
          sprintf(name,"hpulls%d_%d_%d_%d",i,j,k,p);
          hpulls[i][j][k][p] = new TH1F(name,name,50,-5.0,5.0);
        }
      }
    }
  }

  return NOERROR;
}

//------------------------------------------------------------------
// brun
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  //vector<string> factory_names = eventLoop->GetFactoryNames();

  usleep(100000); //this just gives the Main thread a chance to finish printing the "Launching threads" message
  cout<<endl;

  // If int PRINT_ALL is set then add EVERYTHING.
  if(PRINT_ALL){
    //toprint = factory_names;
    SKIP_BORING_EVENTS = 0; // with PRINT_ALL, nothing is boring!
  }else{
    // make sure factories exist for all requested data types
    // If a factory isn't found, but one with a "D" prefixed
    // is, go ahead and correct the name.
    vector<string> really_toprint;
    for(unsigned int i=0; i<toprint.size();i++)
    {
      int found = 0;
      int dfound = 0;
      if(found)
        really_toprint.push_back(toprint[i]);
      else if(dfound)
        really_toprint.push_back("D" + toprint[i]);
      else
        cout<<ansi_red<<"WARNING:"<<ansi_normal
          <<" Couldn't find factory for \""
          <<ansi_bold<<toprint[i]<<ansi_normal
          <<"\"!"<<endl;
    }

    toprint = really_toprint;
  }

  // At this point, toprint should contain a list of all factories
  // in dataClassName:tag format, that both exist and were requested.
  // Seperate the tag from the name and fill the fac_info vector.
  fac_info.clear();
  for(unsigned int i=0;i<toprint.size();i++){
    string name = toprint[i];
    string tag = "";
    unsigned int pos = name.rfind(":",name.size()-1);
    if(pos != (unsigned int)string::npos){
      tag = name.substr(pos+1,name.size());
      name.erase(pos);
    }
    factory_info_t f;
    f.dataClassName = name;
    f.tag = tag;
    fac_info.push_back(f);
  }

  cout<<endl;

  return NOERROR;
}

//------------------------------------------------------------------
// evnt
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{

  if(VERBOSE) cerr << "First thing in evnt....." << endl;

  vector<const DPhoton*> photons;
  vector<DKinematicData> kd_initialState;
  vector<DKinematicData> kd_finalState;
  vector<DKinematicData> kd_initialState_post;
  vector<DKinematicData> kd_finalState_post;
  vector<const DMCThrown*> mcthrown;
  DLorentzVector beam, targ, p, pip[2], pim[2], pi0;
  DLorentzVector beamfit, targfit, pfit, pipfit[2], pimfit[2], pi0fit;
  double mm2, mm2fit, ks[2], ksfit[2];
  double cl;
  double pulls[32];

  DKinematicData kd_beam = DKinematicData();
  DKinematicData kd_targ = DKinematicData();

  // Set the beam
  kd_beam.setMass(0.0);
  kd_beam.setCharge(0);
  kd_beam.setMassFixed();
  kd_beam.setMomentum( TVector3(0.0, 0.0, 9.0) );
  kd_beam.clearErrorMatrix();
  kd_initialState.push_back(kd_beam);
  // Set the target
  kd_targ.setMass(0.93827);
  kd_targ.setCharge(1);
  kd_targ.setMassFixed();
  kd_targ.setMomentum( TVector3(0.0, 0.0, 0.0) );
  kd_targ.clearErrorMatrix();
  kd_initialState.push_back(kd_targ);

  DKinFit *kfit = new DKinFit();
  kfit->SetVerbose(VERBOSE);

  float clCut[5] = {-1.0, 0.01, 0.1, 0.2, 0.5};

  if(eventLoop->Get(mcthrown))
  { 
    kd_finalState.clear();
    int k=0;
    for(int i=0;i<(int)mcthrown.size();i++)
    {
      int pid = mcthrown[i]->type;
      if(VERBOSE) cerr << i << " pid: " << pid << endl;
      if(pid==14||pid==8||pid==9 && i!=3) 
      {
        kd_finalState.push_back(*(mcthrown[i]));
        kd_finalState[k].smearMCThrownMomentum(SMEARWEIGHT);
        k++;
      }
      if(COUNT%100==0) cerr << COUNT << " " << MAX_EVENTS << endl;
      if(COUNT>=MAX_EVENTS) 
      {
        fini();
        exit(1);
      }
    }

    beam = kd_initialState[0].lorentzMomentum();
    targ = kd_initialState[1].lorentzMomentum();
    pip[0] = kd_finalState[0].lorentzMomentum();
    pim[0] = kd_finalState[1].lorentzMomentum();
    pip[1] = kd_finalState[2].lorentzMomentum();
    p = kd_finalState[3].lorentzMomentum();
    pim[1] = beam + targ - p - pip[0] - pim[0] - pip[1];

    mm2 = pim[1].M2();
    ks[0] = (pip[0] + pim[0]).M();
    ks[1] = (pip[1] + pim[1]).M();

    hmm[0][0][0]->Fill(mm2);
    hinvmass[0][0][0][0]->Fill(ks[0]);
    hinvmass[0][0][0][1]->Fill(ks[1]);

    if(VERBOSE) cerr << mm2 << endl;

    // Set up and run the fit
    std::vector<int> constraintParticles;
    constraintParticles.push_back(0); // pi+ 0
    constraintParticles.push_back(1); // pi- 0
    kfit->SetFinal(kd_finalState);
    kfit->SetInitial(kd_initialState);
    kfit->SetMissingParticle(0.13957);
    kfit->SetMassConstraint(0.494,constraintParticles);
    constraintParticles.clear();
    constraintParticles.push_back(2); // pi+ 1
    constraintParticles.push_back(-1); // missing pi- 1
    kfit->SetMassConstraint(0.494,constraintParticles);
    kfit->Fit();

    cl = kfit->Prob();

    hcl[0][0][0]->Fill(cl);

    kd_initialState_post = kfit->GetInitial_out();
    kd_finalState_post = kfit->GetFinal_out();

    beamfit = kd_initialState_post[0].lorentzMomentum();
    targfit = kd_initialState_post[1].lorentzMomentum();
    pipfit[0] = kd_finalState_post[0].lorentzMomentum();
    pimfit[0] = kd_finalState_post[1].lorentzMomentum();
    pipfit[1] = kd_finalState_post[2].lorentzMomentum();
    pfit = kd_finalState_post[3].lorentzMomentum();
    pimfit[1] = beamfit + targfit - pfit - pipfit[0] - pimfit[0] - pipfit[1];

    mm2fit = pimfit[1].M2();
    ksfit[0] = (pipfit[0] + pimfit[0]).M();
    ksfit[1] = (pipfit[1] + pimfit[1]).M();

    for(int k=0;k<12;k++)
    {
      pulls[k] = kfit->GetPull(k+6);
    }
    for(int j=0;j<5;j++)
    {
      if(cl>clCut[j]) 
      {
        hmm[1][0][j]->Fill(mm2fit);
        hinvmass[1][0][j][0]->Fill(ksfit[0]);
        hinvmass[1][0][j][1]->Fill(ksfit[1]);
        for(int k=0;k<12;k++) hpulls[1][0][j][k]->Fill(pulls[k]);
      }
    }

    if(VERBOSE>1)
    {
      for(int k=0;k<2;k++)
      {
        cerr << "initial in: " << kfit->GetInitial_in()[k].energy() << " ";
        cerr << kfit->GetInitial_in()[k].px() << " ";
        cerr << kfit->GetInitial_in()[k].py() << " ";
        cerr << kfit->GetInitial_in()[k].pz() << endl;
        cerr << "initial out: " << kfit->GetInitial_out()[k].energy() << " ";
        cerr << kfit->GetInitial_out()[k].px() << " ";
        cerr << kfit->GetInitial_out()[k].py() << " ";
        cerr << kfit->GetInitial_out()[k].pz() << endl;
      }
      for(int k=0;k<5;k++)
      {
        cerr << "final in: " << kfit->GetFinal_in()[k].energy() << " ";
        cerr << kfit->GetFinal_in()[k].px() << " ";
        cerr << kfit->GetFinal_in()[k].py() << " ";
        cerr << kfit->GetFinal_in()[k].pz() << endl;
        cerr << "final out: " << kfit->GetFinal_out()[k].energy() << " ";
        cerr << kfit->GetFinal_out()[k].px() << " ";
        cerr << kfit->GetFinal_out()[k].py() << " ";
        cerr << kfit->GetFinal_out()[k].pz() << endl;
        int j=k;
        cerr << "mcthrown: " << mcthrown[j]->energy() << " ";
        cerr << mcthrown[j]->px() << " ";
        cerr << mcthrown[j]->py() << " ";
        cerr << mcthrown[j]->pz() << endl;
      }
    }

    if(VERBOSE)
    {
      cerr << "prob after fit: " << kfit->Prob() << endl;
      cerr << "Chi2 after fit: " << kfit->Chi2() << endl;
    }


    COUNT++;
  }

  return NOERROR;
}

//------------------------------------------------------------------
// fini
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
  fout->cd();
  fout->Write();
  fout->cd();
  fout->Close();
  //  delete fout;
  cout<<endl<<"Closed ROOT file"<<endl;

  ///< Called after last event of last event source has been processed.
  return NOERROR;
}				
