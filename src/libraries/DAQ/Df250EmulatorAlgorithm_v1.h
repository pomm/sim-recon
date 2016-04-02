#ifndef _Df250EmulatorAlgorithm_v1_
#define _Df250EmulatorAlgorithm_v1_
#include <JANA/JObject.h>

#include <stdint.h>
#include <vector>
#include <iostream>
using namespace std;

#include <DAQ/Df250EmulatorAlgorithm.h>

/////////////////////////////////////////////////////////////////
// This implements the base class for the f250 firmware emulation
// EmulateFirmware needs to be virtually overwritten by the user
////////////////////////////////////////////////////////////////

class Df250EmulatorAlgorithm_v1:public Df250EmulatorAlgorithm{
    public:

        Df250EmulatorAlgorithm_v1(JEventLoop *loop);
        ~Df250EmulatorAlgorithm_v1(){};

        //Only the emulation routines need to be overwritten
        void EmulateFirmware(const Df250WindowRawData*, Df250PulseTime*, Df250PulsePedestal*, Df250PulseIntegral*);

    protected:
        Df250EmulatorAlgorithm_v1(){};
        // Enables forced use of default values
        int FORCE_DEFAULT;
        // Default values for the essential parameters
        uint32_t NSA_DEF;
        uint32_t NSB_DEF;
        uint16_t THR_DEF;
        // Set verbosity
        int VERBOSE;

};

#endif // _Df250EmulatorAlgorithm_v1_
