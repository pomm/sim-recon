// $Id$
//
//    File: DTOFHit_factory.h
// Created: Wed Aug  7 09:30:17 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTOFHit_factory_
#define _DTOFHit_factory_

#include <vector>
#include <map>
#include <utility>
using namespace std;

#include <JANA/JFactory.h>
#include "DTOFDigiHit.h"
#include "DTOFTDCDigiHit.h"
#include "DTOFHit.h"
#include "DTOFGeometry.h"
#include "TTAB/DTranslationTable.h"
#include "TTAB/DTTabUtilities.h"
using namespace jana;


// store constants so that they can be accessed by plane/bar number
// each entry holds a pair of value for the two different bar ends
// [whatever the conventional choice is]
typedef  vector< vector< pair<double,double> > >  tof_digi_constants_t;


class DTOFHit_factory:public jana::JFactory<DTOFHit>{
 public:
  DTOFHit_factory(){};
  ~DTOFHit_factory(){};
  
  // geometry info - maybe redundant?
  const static int TOF_MAX_CHANNELS = 176;
  int TOF_NUM_PLANES;
  int TOF_NUM_BARS;
  
  // overall scale factors
  double a_scale;
  double t_scale;
  double t_base,t_base_tdc;
  double tdc_adc_time_offset;

  // ADC to Energy conversion for individual PMT channels
  double adc2E[176]; // 4*44 channels
  
  double DELTA_T_ADC_TDC_MAX;
  
  tof_digi_constants_t adc_pedestals;
  tof_digi_constants_t adc_gains;
  tof_digi_constants_t adc_time_offsets;
  tof_digi_constants_t tdc_time_offsets;
  
  vector<vector<double> >timewalk_parameters;
  
  
  DTOFHit* FindMatch(int plane, int bar, int end, double T);
  
  const double GetConstant( const tof_digi_constants_t &the_table,
			    const int in_plane, const int in_bar, 
			    const int in_end ) const;
  const double GetConstant( const tof_digi_constants_t &the_table,
			    const DTOFDigiHit *the_digihit) const;
  const double GetConstant( const tof_digi_constants_t &the_table,
			    const DTOFHit *the_hit) const;
  const double GetConstant( const tof_digi_constants_t &the_table,
			    const DTOFTDCDigiHit *the_digihit) const;
  //const double GetConstant( const tof_digi_constants_t &the_table,
  //			  const DTranslationTable *ttab,
  //			  const int in_rocid, const int in_slot, const int in_channel) const;
  
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  void FillCalibTable(tof_digi_constants_t &table, vector<double> &raw_table,
		      const DTOFGeometry &tofGeom);
  
};

#endif // _DTOFHit_factory_

