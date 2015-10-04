#ifndef __GAUGELINKS__CPP__
#define __GAUGELINKS__CPP__

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              WE DEFINE THE LATTICE GAUGE LINKS AS                                    //
//                          U^{\dagger}_{mu}(x)= P exp[ -ig a_{\mu} A_{\mu}(x)]                         //
//                               WITH DERIVATIVES ACTING AS                                             //
//   \delta U_{mu}(x)/ \delta A_{nu}^{a}(y)= -ig a_{\mu}t^{a} U_{mu}(x) \delta_{x,y} \delta_{\mu\nu}    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class GaugeLinks{
    
public:
    
    INT Volume;
    
    const static INT Dimension=3;
            
    DOUBLE a[Dimension];
    
    DOUBLE aCube;
    
    INT N[Dimension];
        
    SU_Nc_FUNDAMENTAL_FORMAT *U;
    
    INT Index3D(INT x,INT y,INT z){
        
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
        
    }
    
    INT Index(INT x,INT y,INT z,INT mu){
        
        return SUNcGroup::MatrixSize*(mu+Dimension*Index3D(x,y,z));
    }
    
    void GetPosition(INT Index,INT &x,INT &y,INT &z){
        
        z=Index/(N[0]*N[1]); y=(Index-z*(N[0]*N[1]))/(N[0]); x=(Index-N[0]*(y+N[1]*z));
        
    }
    
    DOUBLE* Get(INT x,INT y,INT z,INT mu){
        
        return &(U[Index(x,y,z,mu)]);
    }
    
    void SetIdentity(INT x,INT y,INT z,INT mu){
        
        for(INT z=0;z<=N[2]-1;z++){
            for(INT y=0;y<=N[1]-1;y++){
                for(INT x=0;x<=N[0]-1;x++){
                    for(int mu=0;mu<Dimension;mu++){
                        
                        COPY_SUNcMatrix(this->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                    
                    }
                }
            }
        }
    }

    
    GaugeLinks(INT Nx,INT Ny,INT Nz,INT ax,INT ay,INT az){
        
        this->Volume=Nx*Ny*Nz;
        
        this->N[0]=Nx;
        this->N[1]=Ny;
        this->N[2]=Nz;
        
        this->a[0]=ax;
        this->a[1]=ay;
        this->a[2]=az;
        
        this->aCube=a[0]*a[1]*a[2];
        
        U=new SU_Nc_FUNDAMENTAL_FORMAT[Dimension*SUNcGroup::MatrixSize*Volume];
        
    }
    
    ~GaugeLinks(){
        
        delete U;
    }
    
    
};

void Copy(GaugeLinks *UCopy,GaugeLinks *UOriginal){
    
    for(INT i=0;i<UOriginal->Dimension;i++){
        
        UCopy->a[i]=UOriginal->a[i];
    
        UCopy->N[i]=UOriginal->N[i];
    
    }
    
    UCopy->Volume=UOriginal->Volume;
    
    UCopy->aCube=UOriginal->aCube;
    
    std::memcpy(UCopy->Get(0,0,0,0),UOriginal->Get(0,0,0,0),UOriginal->Dimension*SUNcGroup::MatrixSize*UOriginal->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));

}

void CopyWithBlock(GaugeLinks **UOld,GaugeLinks *UNew){
    
    // COMPUTE BLOCKED GAUGE LINKS //
    for(INT z=0;z<=UNew->N[2]-1;z++){
        for(INT y=0;y<=UNew->N[1]-1;y++){
            for(INT x=0;x<=UNew->N[0]-1;x++){
                
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UNew->Get(x,y,z,0));
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UNew->Get(x,y,z,1));
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UNew->Get(x,y,z,2));
                
            }
        }
    }
    
}

void Block(GaugeLinks **UOld){
    
    // CREATE NEW OBJECTS //
    GaugeLinks *UNew=new GaugeLinks((*UOld)->N[0]/2,(*UOld)->N[1]/2,(*UOld)->N[2]/2,2*(*UOld)->a[0],2*(*UOld)->a[1],2*(*UOld)->a[2]);
    
    // COMPUTE BLOCKED GAUGE LINKS //
    for(INT z=0;z<=UNew->N[2]-1;z++){
        for(INT y=0;y<=UNew->N[1]-1;y++){
            for(INT x=0;x<=UNew->N[0]-1;x++){
                
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UNew->Get(x,y,z,0));
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UNew->Get(x,y,z,1));
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UNew->Get(x,y,z,2));
                
            }
        }
    }
    
    // DELETE OLD OBJECTS //
    delete *UOld;
    
    // SET POINTER TO BLOCKED FIELDS //
    *UOld=UNew;
    
}
void CopyWithRGBlock(GaugeLinks **UOld,GaugeLinks *UNew){
    

    // COMPUTE BLOCKED GAUGE LINKS //
    for(INT z=0;z<=UNew->N[2]-1;z++){
        for(INT y=0;y<=UNew->N[1]-1;y++){
            for(INT x=0;x<=UNew->N[0]-1;x++){
                
                //BUFFERS
                SU_Nc_FUNDAMENTAL_FORMAT UxxP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UyyP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UzzP[SUNcGroup::MatrixSize];
                
                SU_Nc_FUNDAMENTAL_FORMAT UzxP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UzxM[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UyxP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UyxM[SUNcGroup::MatrixSize];
                
                SU_Nc_FUNDAMENTAL_FORMAT UxyP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UxyM[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UzyP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UzyM[SUNcGroup::MatrixSize];
                
                SU_Nc_FUNDAMENTAL_FORMAT UyzP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UyzM[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UxzP[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UxzM[SUNcGroup::MatrixSize];
                
                //SMEARED LINKS FOR BLOCKING
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UxxP);
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UyyP);
                SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UzzP);
                
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,0),(*UOld)->Get(2*x+1,2*y,2*z+1,0),(*UOld)->Get(2*x+2,2*y,2*z,2),UzxP);
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,1),(*UOld)->Get(2*x+1,2*y+1,2*z,1),(*UOld)->Get(2*x,2*y+2,2*z,0),UxyP);
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,2),(*UOld)->Get(2*x,2*y+1,2*z+1,2),(*UOld)->Get(2*x,2*y,2*z+2,1),UyzP);
                
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,0),(*UOld)->Get(2*x+1,2*y+1,2*z,0),(*UOld)->Get(2*x+2,2*y,2*z,1),UyxP);
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,1),(*UOld)->Get(2*x,2*y+1,2*z+1,1),(*UOld)->Get(2*x,2*y+2,2*z,2),UzyP);
                SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,2),(*UOld)->Get(2*x+1,2*y,2*z+1,2),(*UOld)->Get(2*x,2*y,2*z+2,0),UxzP);
                
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y,2*z-1,2),(*UOld)->Get(2*x,2*y,2*z-1,0),(*UOld)->Get(2*x+1,2*y,2*z-1,0),(*UOld)->Get(2*x+2,2*y,2*z-1,2),UzxM);
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x-1,2*y,2*z,0),(*UOld)->Get(2*x-1,2*y,2*z,1),(*UOld)->Get(2*x-1,2*y+1,2*z,1),(*UOld)->Get(2*x-1,2*y+2,2*z,0),UxyM);
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y-1,2*z,1),(*UOld)->Get(2*x,2*y-1,2*z,2),(*UOld)->Get(2*x,2*y-1,2*z+1,2),(*UOld)->Get(2*x,2*y-1,2*z+2,1),UyzM);
                
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y-1,2*z,1),(*UOld)->Get(2*x,2*y-1,2*z,0),(*UOld)->Get(2*x+1,2*y-1,2*z,0),(*UOld)->Get(2*x+2,2*y-1,2*z,1),UyxM);
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y,2*z-1,2),(*UOld)->Get(2*x,2*y,2*z-1,1),(*UOld)->Get(2*x,2*y+1,2*z-1,1),(*UOld)->Get(2*x,2*y+2,2*z-1,2),UzyM);
                SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x-1,2*y,2*z,0),(*UOld)->Get(2*x-1,2*y,2*z,2),(*UOld)->Get(2*x-1,2*y,2*z+1,2),(*UOld)->Get(2*x-1,2*y,2*z+2,0),UxzM);
                
                //RECONSTRUCT BLOCKED LINKS
                SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UxxP,UzxP,UzxM,UyxP,UyxM,UNew->Get(x,y,z,0));
                SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UyyP,UxyP,UxyM,UzyP,UzyM,UNew->Get(x,y,z,1));
                SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UzzP,UyzP,UyzM,UxzP,UxzM,UNew->Get(x,y,z,2));
                
                // RE-UNITARIZE //
                SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,0));
                SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,1));
                SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,2));
                
                
            }
        }
    }
    
}


#endif
