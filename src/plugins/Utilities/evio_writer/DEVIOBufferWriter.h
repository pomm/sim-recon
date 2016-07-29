#ifndef _DEVIOBufferWriter_
#define _DEVIOBufferWriter_

#include <map>
#include <vector>
#include <string>
#include <set>

#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>

#include <JANA/JEventLoop.h>

#include <DAQ/JEventSource_EVIO.h>

#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/Df125WindowRawData.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/Df125Config.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCConfig.h>
#include <DAQ/DCAEN1290TDCConfig.h>
#include <DAQ/DCAEN1290TDCHit.h>
#include <DAQ/DEPICSvalue.h>
#include <DAQ/DEventTag.h>

#include <DANA/DStatusBits.h>
#include <TTAB/DTranslationTable.h>

using namespace std;
using namespace jana;

class DEVIOBufferWriter
{
  public:
    DEVIOBufferWriter(bool compact_flag = false, bool prefer_emulated_flag = false) {
        COMPACT = compact_flag;
        PREFER_EMULATED = prefer_emulated_flag;

        write_out_all_rocs = true;   // default to writing all data
    }
    ~DEVIOBufferWriter(void);

    void WriteEventToBuffer(JEventLoop *loop, vector<uint32_t> &buff, vector<const JObject *> objects_to_save) const;
    void WriteEventToBuffer(JEventLoop *locEventLoop, vector<uint32_t> &buff) const;

    void SetROCsToWriteOut(set<uint32_t> &new_rocs_to_write_out) {
        rocs_to_write_out = new_rocs_to_write_out;

        if(rocs_to_write_out.size() == 0)
            write_out_all_rocs = true;
    }


  protected:

        void WriteBuiltTriggerBank(vector<uint32_t> &buff,
                                   JEventLoop *loop,
                                   vector<const DCODAROCInfo*> &coda_rocinfos,
                                   vector<const DCODAEventInfo*> &coda_events) const;

		void WriteCAEN1290Data(vector<uint32_t> &buff,
                               vector<const DCAEN1290TDCHit*>    &caen1290hits,
                               vector<const DCAEN1290TDCConfig*> &caen1290configs,
                               unsigned int Nevents) const;


		void WriteF1Data(vector<uint32_t> &buff,
                         vector<const DF1TDCHit*>          &F1hits,
                         vector<const DF1TDCTriggerTime*>  &F1tts,
                         vector<const DF1TDCConfig*>       &F1configs,
                         unsigned int Nevents) const;

		void Writef250Data(vector<uint32_t> &buff,
                           vector<const Df250PulseIntegral*> &f250pis,
                           vector<const Df250TriggerTime*>   &f250tts,
                           vector<const Df250WindowRawData*> &f250wrds,
                           unsigned int Nevents) const;

		void Writef125Data(vector<uint32_t> &buff,
                           vector<const Df125PulseIntegral*> &f125pis,
                           vector<const Df125CDCPulse*>      &f125cdcpulses,
                           vector<const Df125FDCPulse*>      &f125fdcpulses,
                           vector<const Df125TriggerTime*>   &f125tts,
                           vector<const Df125WindowRawData*> &f125wrds,
                           vector<const Df125Config*>        &f125configs,
                           unsigned int Nevents) const;

		void WriteEPICSData(vector<uint32_t> &buff,
                            vector<const DEPICSvalue*> epicsValues) const;

		void WriteEventTagData(vector<uint32_t> &buff,
                               uint64_t event_status,
                               const DL3Trigger* l3trigger) const;

        void WriteBORData(JEventLoop *loop,
                          vector<uint32_t> &buff) const;

        void WriteTSSyncData(JEventLoop *loop,
                             vector<uint32_t> &buff,
                             const DL1Info *l1info) const;

        bool write_out_all_rocs;
        set<uint32_t> rocs_to_write_out;

        bool COMPACT;
        bool PREFER_EMULATED;
};

#endif
