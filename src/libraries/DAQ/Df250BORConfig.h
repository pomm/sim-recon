// $Id$
//
//    File: Df250BORConfig.h
// Created: Tue Jan 26 13:04:46 EST 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _Df250BORConfig_
#define _Df250BORConfig_

#include <JANA/JObject.h>

#include <DAQ/bor_roc.h>

// This class inherits both from JObject and f250config. The former
// so that it can be incorporated easily into the JANA framework.
// The latter so we can use the data struct defined in bor_roc.h.
// The file bor_roc.h exists in 2 places:
//
//  1. in the DAQ library of sim-recon
//  2. in the vme/src/rcm/monitor directory in the online
//


class Df250BORConfig:public jana::JObject, public f250config{
	public:
		JOBJECT_PUBLIC(Df250BORConfig);

		Df250BORConfig(){}
		virtual ~Df250BORConfig(){}
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid"     , "%d", rocid);
			AddString(items, "slot"      , "%d", slot);
			AddString(items, "version"   , "0x%x", version);
			AddString(items, "ctrl1"     , "0x%x", ctrl1);
			AddString(items, "ctrl2"     , "0x%x", ctrl2);
			AddString(items, "blk_level" , "%d", blk_level);
			AddString(items, "ptw"       , "%d", adc_ptw);
			AddString(items, "pl"        , "%d", adc_pl);
			AddString(items, "nsb"       , "%d", adc_nsb);
			AddString(items, "nsa"       , "%d", adc_nsa);
		}

};

#endif // _Df250BORConfig_

