#ifndef __LATTICE__CPP__
#define __LATTICE__CPP__

///////////////////////////////////////////////
//CONTAINS LATTICE DISCRETIZATION PARAMETERS //
///////////////////////////////////////////////

namespace Lattice{
    
    //DIMENSION OF THE LATTICE
    static const int Dimension=3;
    
    //NUMBER OF INITIAL LATTICE SITES AND LATTICE SPACINGS
    INT N[3]={64,64,64};
    static const DOUBLE a[3]={1.0,1.0,1.0};
    
    // SCALE //
    static DOUBLE aScale=1.0;
    
    //LATTICE VOLUME
    INT Volume=64*64*64;
    
}
    
#endif
