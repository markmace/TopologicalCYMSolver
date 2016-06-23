#ifndef __GAUGETRANSFORMATIONS__CPP__
#define __GAUGETRANSFORMATIONS__CPP__

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WE DEFINE THE LATTICE GAUGE LINKS AS                                    //
//                          U^{\dagger}_{mu}(x)= P exp[ -ig a_{\mu} A_{\mu}(x)]                         //
//                               WITH DERIVATIVES ACTING AS                                             //
//   \delta U_{mu}(x)/ \delta A_{nu}^{a}(y)= -ig a_{\mu}t^{a} U_{mu}(x) \delta_{x,y} \delta_{\mu\nu}    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class GaugeTransformations{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    DOUBLE a[Dimension];
    
    INT N[Dimension];
        
    SU_Nc_FUNDAMENTAL_FORMAT *G;
    
    INT Index3D(INT x,INT y,INT z){
        
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
        
    }
    
    INT Index(INT x,INT y,INT z){

        return SUNcGroup::MatrixSize*(Index3D(x,y,z));
    
    }
    
    SU_Nc_FUNDAMENTAL_FORMAT* Get(INT x,INT y,INT z){
        
        return &(G[Index(x,y,z)]);
    }
    
    GaugeTransformations(INT Nx,INT Ny,INT Nz,INT ax,INT ay,INT az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        G=new SU_Nc_FUNDAMENTAL_FORMAT[SUNcGroup::MatrixSize*Volume];
        
    }
    
    ~GaugeTransformations(){
        
        delete G;
    }
    
    
};

void Save(std::string fname,GaugeTransformations *GaugeTrafo){
    
    std::ofstream OutStream;
    
    OutStream.open(fname.c_str());
    
    for(INT z=0;z<GaugeTrafo->N[2];z++){
        for(INT y=0;y<GaugeTrafo->N[1];y++){
            for(INT x=0;x<GaugeTrafo->N[0];x++){
                
                OutStream << x << " " << y << " " << z << " " << SUNcGroup::IO::MatrixToString(GaugeTrafo->Get(x,y,z)) << std::endl;
                
            }
        }
    }
    
    OutStream.close();
    
    
}

#endif