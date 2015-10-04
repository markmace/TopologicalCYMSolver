#ifndef __GLINKSEFIELDS__CPP__
#define __GLINKSEFIELDS__CPP__

namespace GLinks{
    
    GaugeLinks *U;
    
    void Init(){
        
        U=new GaugeLinks(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);

    }
    
}

namespace EFields{
    
    ElectricFields *E;
    
    void Init(){
        
        E=new ElectricFields(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);

    }
    
}

namespace Lattice{
    
    ///////////////////////////
    //INITIALIZATION ROUTINE //
    ///////////////////////////
   
    DOUBLE aCube;
    
    void Init(){
        
        GLinks::Init();
        
        EFields::Init();
        
        aCube=GLinks::U->a[0]*GLinks::U->a[1]*GLinks::U->a[2];

        
    }
    

}
#endif