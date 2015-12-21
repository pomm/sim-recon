#if !defined(TWOPISODING)
#define TWOPISODING

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPUTwoPiSoding_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                          GDouble mass0, GDouble width0, int orbitL,
                          int daught1, int daught2 );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class TwoPiSoding : public UserAmplitude< TwoPiSoding >
{
  
public:
	
	TwoPiSoding() : UserAmplitude< TwoPiSoding >() {}
	TwoPiSoding( const vector< string >& args );
	
  ~TwoPiSoding(){}
  
	string name() const { return "TwoPiSoding"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_mass0;
  AmpParameter m_width0;
  int m_orbitL;
  
  pair< string, string > m_daughters;  

  AmpParameter m_soding;
  AmpParameter m_slope;

  AmpParameter rho000;
  AmpParameter rho100;
  AmpParameter rho1m10;

  AmpParameter rho111;
  AmpParameter rho001;
  AmpParameter rho101;
  AmpParameter rho1m11;

  AmpParameter rho102;
  AmpParameter rho1m12;
};

#endif
