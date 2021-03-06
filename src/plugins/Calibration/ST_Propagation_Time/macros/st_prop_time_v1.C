// File: st_tw_fits.C
// Last Modified: 01/12/2016
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms and propagation time correction.
#include "TH1.h"
#include <TH2I.h>
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <TDirectory.h>
#include <TLine.h>
#include <TPaveLabel.h>
#include <TPaletteAxis.h>
#include <stdio.h>
#include <stdint.h>
#include <fstream>
using namespace std;
using std::string;
// ***************** define constants and varibles*************************
const Int_t NCHANNELS          = 30;
Double_t t_z_ss_fit_slope[NCHANNELS][3];
Double_t t_z_ss_fit_slope_err[NCHANNELS][3];
Double_t t_z_ss_fit_intercept[NCHANNELS][3];
Double_t t_z_ss_fit_intercept_err[NCHANNELS][3];

Double_t t_z_bs_fit_slope[NCHANNELS][3];
Double_t t_z_bs_fit_slope_err[NCHANNELS][3];
Double_t t_z_bs_fit_intercept[NCHANNELS][3];
Double_t t_z_bs_fit_intercept_err[NCHANNELS][3];

Double_t t_z_ns_fit_slope[NCHANNELS][3];
Double_t t_z_ns_fit_slope_err[NCHANNELS][3];
Double_t t_z_ns_fit_intercept[NCHANNELS][3];
Double_t t_z_ns_fit_intercept_err[NCHANNELS][3];
const int NOXbins = 1300;
const int NOYbins = 200;
int sss[NCHANNELS][NOXbins][NOYbins];
Double_t binsss[NCHANNELS][NOXbins][NOYbins];
Double_t ss_max[NCHANNELS];
Double_t low_cut = 0.25;
// Declare canvas
TCanvas *PT_can[30];

void st_prop_time_v1(char*input_filename)
//void st_tw_fits()
{
  TFile *df = new TFile(input_filename);
  std::ofstream st_prop_time_constants;
  st_prop_time_constants.open ("st_prop_time_constants.txt", std::ofstream::out);
  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
  for (unsigned int j = 0; j < NCHANNELS; j++)
    {
      //Create the canvas
      PT_can[j] = new TCanvas( Form("PT_can_%i",j+1),  Form("PT_can_%i",j+1), 800, 450);
      PT_can[j]->Divide(3, 1);
	  
      // Grab the histograms 
      char* ss = Form("h2_PropTime_z_SS_chan_%i",j+1);
      TH1I* h2_ss = (TH1I*) df->Get(ss);
      char* bs = Form("h2_PropTime_z_BS_chan_%i",j+1);
      TH1I* h2_bs = (TH1I*) df->Get(bs);
      char* ns = Form("h2_PropTime_z_NS_chan_%i",j+1);
      TH1I* h2_ns = (TH1I*) df->Get(ns);
      
      cout << "==================================================" << endl;
      cout << "Processing Channel " << j+1 << endl;
      cout << "==================================================" << endl;
      //*****************************************************************
      //***********  Plot stt vs z SS                  ******************
      //*****************************************************************
      PT_can[j]->cd(1);
      gStyle->SetOptFit(111111);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h2_ss->SetTitle("SS ST Time Vs Z");
      h2_ss->GetZaxis()->SetRangeUser(0,4);
      h2_ss->GetYaxis()->SetRangeUser(0.,4.);
      h2_ss->GetYaxis()->SetTitle("ST Time (ns)");
      h2_ss->GetXaxis()->SetRangeUser(50.,78.0);
      h2_ss->Draw("colz");
      h2_ss->Fit("pol1");		
      t_z_ss_fit_slope[j][1] = pol1->GetParameter(1);
      t_z_ss_fit_slope_err[j][1] = pol1->GetParError(1);
      t_z_ss_fit_intercept[j][0] = pol1->GetParameter(0);
      t_z_ss_fit_intercept_err[0][1] = pol1->GetParError(0);

      //*****************************************************************
      //***********  Plot stt vs z BS                  ******************
      //*****************************************************************
      PT_can[j]->cd(2);
      gPad->SetTicks();
      gPad->SetGrid();
      h2_bs->Draw("colz"); 
      h2_bs->SetTitle("BS ST Time Vs Z");
      h2_bs->GetZaxis()->SetRangeUser(0,4);
      h2_bs->GetYaxis()->SetRangeUser(0.,4.);
      h2_bs->GetYaxis()->SetTitle("ST Time (ns)");
      h2_bs->GetXaxis()->SetRangeUser(77.,81.0);
     

      h2_bs->Fit("pol1");
      t_z_bs_fit_slope[j][1] = pol1->GetParameter(1);
      t_z_bs_fit_slope_err[j][1] = pol1->GetParError(1);
      t_z_bs_fit_intercept[j][0] = pol1->GetParameter(0);
      t_z_bs_fit_intercept_err[0][1] = pol1->GetParError(0);
      //*****************************************************************
      //***********  Plot stt vs z NS                  ******************
      //*****************************************************************
      PT_can[j]->cd(3);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h2_ns->Draw("colz");
      h2_ns->SetTitle("NS ST Time Vs Z");
      h2_ns->GetZaxis()->SetRangeUser(0,4);
      h2_ns->GetYaxis()->SetRangeUser(0.,4.);
      h2_ns->GetYaxis()->SetTitle("ST Time (ns)");
      h2_ns->GetXaxis()->SetRangeUser(80.,98.0);
      h2_ns->Fit("pol1");
      t_z_ns_fit_slope[j][1] = pol1->GetParameter(1);
      t_z_ns_fit_slope_err[j][1] = pol1->GetParError(1);
      t_z_ns_fit_intercept[j][0] = pol1->GetParameter(0);
      t_z_ns_fit_intercept_err[0][1] = pol1->GetParError(0);
      //*****************************************************************
      //***********  save data to dump it to ccdb         ******************
      //*****************************************************************
      st_prop_time_constants << "\t"<< t_z_ss_fit_intercept[j][0] << "\t"<< t_z_ss_fit_slope[j][1]<< "\t"<< t_z_bs_fit_intercept[j][0] << "\t"<< t_z_bs_fit_slope[j][1]<< "\t"<< t_z_ns_fit_intercept[j][0] << "\t"<< t_z_ns_fit_slope[j][1]<< endl;     
    }
  st_prop_time_constants.close();
}

