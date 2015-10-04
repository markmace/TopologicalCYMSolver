#ifndef __FOURIER_SPACE_CPP__
#define __FOURIER_SPACE_CPP__

// INCLUDE THREE DIMENSIONAL FAST FOURIER TRANSFORM //
#include "../FFT/FFT3D.cpp"

namespace FourierSpace{
        
    FFT3D *E0;
    FFT3D *E1;
    FFT3D *E2;
    
    FFT3D *A0;
    FFT3D *A1;
    FFT3D *A2;
    
    FFT3D *GaugeUpdate;
    FFT3D *BlockedGaugeUpdate;

    // INITIALIZATION //
    void Init(){
        
        A0=new FFT3D(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],SUNcAlgebra::VectorSize);
        A1=new FFT3D(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],SUNcAlgebra::VectorSize);
        A2=new FFT3D(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],SUNcAlgebra::VectorSize);
        
        E0=new FFT3D(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],SUNcAlgebra::VectorSize);
        E1=new FFT3D(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],SUNcAlgebra::VectorSize);
        E2=new FFT3D(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],SUNcAlgebra::VectorSize);
        
        
        GaugeUpdate=new FFT3D(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],SUNcAlgebra::VectorSize);
        BlockedGaugeUpdate=new FFT3D(GLinks::U->N[0]/2,GLinks::U->N[1]/2,GLinks::U->N[2]/2,SUNcAlgebra::VectorSize);

        
    }
    
}

#endif
