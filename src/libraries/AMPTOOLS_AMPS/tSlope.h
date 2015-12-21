#if !defined(TSLOPE)
#define TSLOPE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUtSlope_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     int j, int m, GDouble bigTheta, GDouble refFact );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class tSlope : public UserAmplitude< tSlope >
{
    
public:
	
	tSlope() : UserAmplitude< tSlope >() { };
	tSlope( const vector< string >& args );
	
	string name() const { return "tSlope"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:

  AmpParameter slope;

  vector<int> resonanceChildren;

};

#endif
