#if !defined(TWOPIRHOJPAC)
#define TWOPIRHOJPAC

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUTwoPiRhoJPAC_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     int j, int m, GDouble bigTheta, GDouble refFact );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

struct dcomplex {
	double dr;
	double di;
};

class TwoPiRhoJPAC : public UserAmplitude< TwoPiRhoJPAC >
{
    
public:
	
	TwoPiRhoJPAC() : UserAmplitude< TwoPiRhoJPAC >() { };
	TwoPiRhoJPAC( const vector< string >& args );
	
	string name() const { return "TwoPiRhoJPAC"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;

#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:

  AmpParameter param0;
  AmpParameter param1;
  AmpParameter param2;
  AmpParameter param3;

  AmpParameter b0;
  AmpParameter b1;
  AmpParameter b2;
  AmpParameter b3;

};

#endif
