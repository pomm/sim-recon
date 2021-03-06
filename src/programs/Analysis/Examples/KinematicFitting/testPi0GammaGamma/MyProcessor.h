// $Id: MyProcessor.h 2319 2006-12-12 18:23:09Z davidl $
// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// hd_dump print event info to screen
///

#ifndef _MYPROCESSOR_
#define _MYPROCESSOR_


#include "JANA/JEventProcessor.h"
#include "JANA/JEventLoop.h"
#include "JANA/JFactory.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

extern int PAUSE_BETWEEN_EVENTS;
extern int SKIP_BORING_EVENTS;
extern int PRINT_ALL;
extern float MAX_EVENTS;
extern int COUNT;
extern int VERBOSE;
extern char* OUTNAME;
extern float MASS;
extern float ERRMATRIXWEIGHT;
//extern TFile *fout;
//extern TH1F *hpi0[4];
//extern TTree *tree;

extern vector<string> toprint;

class MyProcessor:public JEventProcessor
{
	public:
		jerror_t init(void);				///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventnumber);						///< Called every event.
		jerror_t erun(void){return NOERROR;};				///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);				///< Called after last event of last event source has been processed.

		typedef struct{
			string dataClassName;
			string tag;
		}factory_info_t;
		vector<factory_info_t> fac_info;

  private:
    TFile* fout;
    TH1F* hpi0[4][4][5];
    TH1F* hpi0fit[4][4][5];
    TH1F* hcl[4][4][5];
    TH1F* hclGood[4];
    TH1F* hpulls[4][4][5][6];
};

#endif
