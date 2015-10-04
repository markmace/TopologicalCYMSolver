#ifndef __ELECTRICFIELDS__CPP__
#define __ELECTRICFIELDS__CPP__

////////////////////////////////////////////////////////////////////////////////////////////////
//                    WE DEFINE THE DIMENSIONLESS ELECTRIC FIELD VARIABLES AS                 //
//                 E^{\mu}(x)=g \sqrt{-g(x)} a^3/a_{mu}  g^{\mu\nu} \partial_{\tau} A_{\nu}   //
////////////////////////////////////////////////////////////////////////////////////////////////

class ElectricFields{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
    
    // FOR SU(2), SHOULD MAKE FLAG //
    //const static INT VectorSize=3;
    
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];
        
    SU_Nc_ALGEBRA_FORMAT *E;
    
    // FOR 3D //
    
    INT Index3D(INT x,INT y,INT z){
        
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
        
    }
    
    INT Index(INT x,INT y,INT z,INT mu,INT a){
        
        //((a)+SUNcAlgebra::VectorSize*((mu)+Lattice::Dimension*Index3D((x),(y),(z))))
        return a+SUNcAlgebra::VectorSize*(mu+Dimension*Index3D(x,y,z));
        
    }
    
    DOUBLE* Get(INT x,INT y,INT z,INT mu,INT a){
        
        return &E[Index(x,y,z,mu,a)];
    }

    void SetZero(){
        
        for(INT z=0;z<=N[2]-1;z++){
            for(INT y=0;y<=N[1]-1;y++){
                for(INT x=0;x<=N[0]-1;x++){
                    for(int mu=0;mu<Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            this->Get(x,y,z,mu,a)[0]=DOUBLE(0.0);
                        
                        }
                    }
                }
            }
        }
    }
    
    
    ElectricFields(INT Nx,INT Ny,INT Nz,INT ax,INT ay,INT az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        this->aCube=a[0]*a[1]*a[2];
        
        E=new SU_Nc_ALGEBRA_FORMAT[Dimension*SUNcAlgebra::VectorSize*Volume];
        
    }
    
    ~ElectricFields(){
        
        delete E;
    }
};

void Copy(ElectricFields *ECopy,ElectricFields *EOriginal){
    
    for(INT i=0;i<EOriginal->Dimension;i++){
        
        ECopy->a[i]=EOriginal->a[i];
        
        ECopy->N[i]=EOriginal->N[i];
        
    }
    
    ECopy->Volume=EOriginal->Volume;
    
    ECopy->aCube=EOriginal->aCube;
    
    std::memcpy(ECopy->Get(0,0,0,0,0),EOriginal->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::VectorSize*Lattice::Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
    
    
}

/////////////////////////
//BLOCK ELECTRIC FIELDS//
/////////////////////////

void Block(ElectricFields **EOld){
    
    // CREATE NEW OBJECT //
    ElectricFields *ENew=new ElectricFields((*EOld)->N[0]/2,(*EOld)->N[1]/2,(*EOld)->N[2]/2,2*(*EOld)->a[0],2*(*EOld)->a[1],2*(*EOld)->a[2]);
    
    // SET ELECTRIC FIELDS TO ZERO //
    ENew->SetZero();
    
    // DELETE OLD OBJECT //
    delete *EOld;
    
    // SET POINTER TO BLOCKED FIELDS //
    *EOld=ENew;
    
}

#endif

