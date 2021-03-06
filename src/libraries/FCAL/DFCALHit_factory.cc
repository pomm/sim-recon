// $Id$
//
//    File: DFCALHit_factory.cc
// Created: Tue Aug  6 12:23:43 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALHit_factory.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250PulsePedestal.h"
#include "DAQ/Df250Config.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DFCALHit_factory::init(void)
{
    // initialize calibration tables
    vector< vector<double > > new_gains(DFCALGeometry::kBlocksTall, 
            vector<double>(DFCALGeometry::kBlocksWide));
    vector< vector<double > > new_pedestals(DFCALGeometry::kBlocksTall, 
            vector<double>(DFCALGeometry::kBlocksWide));
    vector< vector<double > > new_t0s(DFCALGeometry::kBlocksTall, 
            vector<double>(DFCALGeometry::kBlocksWide));
    vector< vector<double > > new_qualities(DFCALGeometry::kBlocksTall, 
            vector<double>(DFCALGeometry::kBlocksWide));

    gains = new_gains;
    pedestals = new_pedestals;
    time_offsets = new_t0s;
    block_qualities = new_qualities;

    // set the base conversion scales --
    // a_scale should definitely come from
    // the DB so set it to a value that will
    // be noticeably wrong
    a_scale = 0.0;      // GeV/FADC integral unit
    t_scale = 0.0625;   // 62.5 ps/count
    t_base  = 0.;       // ns

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFCALHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    // Only print messages for one thread whenever run number change
    static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    static set<int> runs_announced;
    pthread_mutex_lock(&print_mutex);
    bool print_messages = false;
    if(runs_announced.find(runnumber) == runs_announced.end()){
        print_messages = true;
        runs_announced.insert(runnumber);
    }
    pthread_mutex_unlock(&print_mutex);

    // extract the FCAL Geometry
    vector<const DFCALGeometry*> fcalGeomVect;
    eventLoop->Get( fcalGeomVect );
    if (fcalGeomVect.size() < 1)
        return OBJECT_NOT_AVAILABLE;
    const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

    /// Read in calibration constants
    vector< double > raw_gains;
    vector< double > raw_pedestals;
    vector< double > raw_time_offsets;
    vector< double > raw_block_qualities;    // we should change this to an int?

    if(print_messages) jout << "In DFCALHit_factory, loading constants..." << endl;

    // load scale factors
    map<string,double> scale_factors;
    if (eventLoop->GetCalib("/FCAL/digi_scales", scale_factors))
        jout << "Error loading /FCAL/digi_scales !" << endl;
    if (scale_factors.find("FCAL_ADC_ASCALE") != scale_factors.end())
        a_scale = scale_factors["FCAL_ADC_ASCALE"];
    else
        jerr << "Unable to get FCAL_ADC_ASCALE from /FCAL/digi_scales !" << endl;
    if (scale_factors.find("FCAL_ADC_TSCALE") != scale_factors.end())
        t_scale = scale_factors["FCAL_ADC_TSCALE"];
    else
        jerr << "Unable to get FCAL_ADC_TSCALE from /FCAL/digi_scales !" << endl;

    // load base time offset
    map<string,double> base_time_offset;
    if (eventLoop->GetCalib("/FCAL/base_time_offset",base_time_offset))
        jout << "Error loading /FCAL/base_time_offset !" << endl;
    if (base_time_offset.find("FCAL_BASE_TIME_OFFSET") != base_time_offset.end())
        t_base = base_time_offset["FCAL_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get FCAL_BASE_TIME_OFFSET from /FCAL/base_time_offset !" << endl;

    // load constant tables
    if (eventLoop->GetCalib("/FCAL/gains", raw_gains))
        jout << "Error loading /FCAL/gains !" << endl;
    if (eventLoop->GetCalib("/FCAL/pedestals", raw_pedestals))
        jout << "Error loading /FCAL/pedestals !" << endl;
    if (eventLoop->GetCalib("/FCAL/timing_offsets", raw_time_offsets))
        jout << "Error loading /FCAL/timing_offsets !" << endl;
    if (eventLoop->GetCalib("/FCAL/block_quality", raw_block_qualities))
        jout << "Error loading /FCAL/block_quality !" << endl;

    FillCalibTable(gains, raw_gains, fcalGeom);
    FillCalibTable(pedestals, raw_pedestals, fcalGeom);
    FillCalibTable(time_offsets, raw_time_offsets, fcalGeom);
    FillCalibTable(block_qualities, raw_block_qualities, fcalGeom);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DFCALHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    /// Generate DFCALHit object for each DFCALDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.
    char str[256];

    // extract the FCAL Geometry (for positionOnFace())
    vector<const DFCALGeometry*> fcalGeomVect;
    eventLoop->Get( fcalGeomVect );
    if (fcalGeomVect.size() < 1)
        return OBJECT_NOT_AVAILABLE;
    const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

    vector<const DFCALDigiHit*> digihits;
    loop->Get(digihits);
    for (unsigned int i=0; i < digihits.size(); i++) {

        const DFCALDigiHit *digihit = digihits[i];

        // There is a slight difference between Mode 7 and 8 data
        // The following condition signals an error state in the flash algorithm in both modes
        // Do not make hits out of these
        const Df250PulsePedestal* PPobj = NULL;
        digihit->GetSingle(PPobj);
        if (PPobj != NULL){
            if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
        }
        //if (digihit->pulse_time==0) continue;

        // Check to see if the hit corresponds to a valid channel
        if (fcalGeom.isBlockActive(digihit->row,digihit->column) == false) {
            sprintf(str, "DFCALHit corresponds to inactive channel!  "
                    "row=%d, col=%d", 
                    digihit->row, digihit->column);
            throw JException(str);
        }

        // throw away hits from bad or noisy channels
        fcal_quality_state quality = 
            static_cast<fcal_quality_state>(block_qualities[digihit->row][digihit->column]);
        if ( (quality==BAD) || (quality==NOISY) ) continue;

        // get pedestal from CCDB -- we should use it instead
        // of the event-by-even pedestal

        double pedestal = pedestals[digihit->row][digihit->column];
        double integratedPedestal = 0.0;
        const Df250PulseIntegral* PIobj = NULL;

        digihit->GetSingle(PIobj);

        if( PIobj != NULL ){

            if( pedestal == 0 ) {

                // we should use the fixed database pedestal
                // object as it is less susceptible to noise
                // than the event-by-event pedestal

                // if the database pedestal is zero then try
                // the event-by-event one:

                pedestal = (double)PIobj->pedestal / 
                    (double)PIobj->nsamples_pedestal;
            }

            double nsamples_integral = (double)PIobj->nsamples_integral;
            integratedPedestal = pedestal * nsamples_integral;
        }
        else{

	  static uint64_t nWarningsFCALIntegral = 0;
	  if(++nWarningsFCALIntegral <= 20) cerr << "ERROR! no associated FCAL integral object." << endl;
	  if(nWarningsFCALIntegral == 20  ) cerr << "DFCALHit_factory: LAST WARNING (others will be suppressed)" << endl;
        }

        double pulse_amplitude = 0;

        if( PPobj != NULL ){

            pulse_amplitude = (double)PPobj->pulse_peak - pedestal;
        }
        else{

	  static uint64_t nWarningsFCALPedestal = 0;
	  if(++nWarningsFCALPedestal <= 20) cerr << "ERROR! no associated FCAL pedestal object." << endl;
	  if(nWarningsFCALPedestal == 20  ) cerr << "DFCALHit_factory: LAST WARNING (others will be suppressed)" << endl;
	  
        }

        DFCALHit *hit = new DFCALHit;
        hit->row    = digihit->row;
        hit->column = digihit->column;

        // Apply calibration constants
        double A = (double)digihit->pulse_integral;
        double T = (double)digihit->pulse_time;
        hit->E = a_scale * gains[hit->row][hit->column] * (A - integratedPedestal);
        hit->t = t_scale * T - time_offsets[hit->row][hit->column] + t_base;

        // Get position of blocks on front face. (This should really come from
        // hdgeant directly so the poisitions can be shifted in mcsmear.)
        DVector2 pos = fcalGeom.positionOnFace(hit->row, hit->column);
        hit->x = pos.X();
        hit->y = pos.Y();

        // recored the pulse integral to peak ratio since this is
        // a useful quality metric for the PMT pulse
        hit->intOverPeak = ( A - integratedPedestal ) / pulse_amplitude;

        // do some basic quality checks before creating the objects
        if( ( hit->E > 0 ) &&
                ( digihit->pulse_time > 0 ) 
          ){

            hit->AddAssociatedObject(digihit);
            _data.push_back(hit);
        }
        else{

            delete hit;
        }
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DFCALHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DFCALHit_factory::fini(void)
{
    return NOERROR;
}

//------------------
// FillCalibTable
//------------------
void DFCALHit_factory::FillCalibTable( fcal_digi_constants_t &table, 
        const vector<double> &raw_table, 
        const DFCALGeometry &fcalGeom)
{
    char str[256];

    // sanity check that we have the right geometry
    // (deprecate this?) 
    if (fcalGeom.numActiveBlocks() != FCAL_MAX_CHANNELS) {
        sprintf(str, "FCAL geometry is wrong size! channels=%d (should be %d)", 
                fcalGeom.numActiveBlocks(), FCAL_MAX_CHANNELS);
        throw JException(str);
    }

    // check to see if the table is the right size
    if ( fcalGeom.numActiveBlocks() != static_cast<int>(raw_table.size()) ) {
        sprintf(str, "FCAL constant table is wrong size! channels=%d (should be %d)", 
                fcalGeom.numActiveBlocks(), static_cast<int>(raw_table.size()));
        throw JException(str);
    }

    for (int channel=0; channel < static_cast<int>(raw_table.size()); channel++)
    {
        // make sure that we don't try to load info for channels that don't exist
        if (channel == fcalGeom.numActiveBlocks())
            break;

        int row = fcalGeom.row(channel);
        int col = fcalGeom.column(channel);

        // results from DFCALGeometry should be self consistent, but add in some
        // sanity checking just to be sure
        if (fcalGeom.isBlockActive(row,col) == false) {
            sprintf(str, "Loading FCAL constant for inactive channel!  "
                    "row=%d, col=%d", row, col);
            throw JException(str);
        }

        table[row][col] = raw_table[channel];
    }
}

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
        const int in_row, 
        const int in_column) const
{
    char str[256];

    if ( (in_row <= 0) || (in_row > DFCALGeometry::kBlocksTall)) {
        sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_row, DFCALGeometry::kBlocksTall);
        cerr << str << endl;
        throw JException(str);
    }
    if ( (in_column <= 0) || (in_column > DFCALGeometry::kBlocksWide)) {
        sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_column, DFCALGeometry::kBlocksWide);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_row][in_column];
}

const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
        const DFCALDigiHit *in_digihit) const
{
    char str[256];

    if ( (in_digihit->row <= 0) || (in_digihit->row > DFCALGeometry::kBlocksTall)) {
        sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()!"
                " requested=%d , should be %ud", 
                in_digihit->row, DFCALGeometry::kBlocksTall);
        cerr << str << endl;
        throw JException(str);
    }
    if ( (in_digihit->column <= 0) || (in_digihit->column > DFCALGeometry::kBlocksWide)) {
        sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()!"
                " requested=%d , should be %ud",
                in_digihit->column, DFCALGeometry::kBlocksWide);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_digihit->row][in_digihit->column];
}

const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
        const DFCALHit *in_hit) const {

    char str[256];

    if ( (in_hit->row <= 0) || (in_hit->row > DFCALGeometry::kBlocksTall)) {
        sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()! "
                "requested=%d , should be %ud", in_hit->row, DFCALGeometry::kBlocksTall);
        cerr << str << endl;
        throw JException(str);
    }
    if ( (in_hit->column <= 0) || (in_hit->column > DFCALGeometry::kBlocksWide)) {
        sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()!"
                " requested=%d , should be %ud", 
                in_hit->column, DFCALGeometry::kBlocksWide);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_hit->row][in_hit->column];
}
/*
   const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
   const DTranslationTable *ttab,
   const int in_rocid,
   const int in_slot,
   const int in_channel) const
   {
   char str[256];

   DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
   DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);

   if ( (channel_info.fcal.row <= 0) 
   || (channel_info.fcal.row > static_cast<unsigned int>(DFCALGeometry::kBlocksTall))) {
   sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()!"
   " requested=%d , should be %ud", 
   channel_info.fcal.row, DFCALGeometry::kBlocksTall);
   cerr << str << endl;
   throw JException(str);
   }
   if ( (channel_info.fcal.col <= 0) 
   || (channel_info.fcal.col > static_cast<unsigned int>(DFCALGeometry::kBlocksWide))) {
   sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()!"
   " requested=%d , should be %ud",
   channel_info.fcal.row, DFCALGeometry::kBlocksWide);
   cerr << str << endl;
   throw JException(str);
   }

   return the_table[channel_info.fcal.row][channel_info.fcal.col];
   }
   */
