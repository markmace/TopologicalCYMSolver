#ifndef __GAUGETRANSFORMATION__CPP__
#define __GAUGETRANSFORMATION__CPP__

////////////////////////////////////////////////
//GAUGE TRANSFORMED LINKS AND ELECTRIC FIELDS //
////////////////////////////////////////////////

namespace GaugeFixedVariables{
    
    GaugeLinks *U;
    ElectricFields *E;
    
    //INITIALIZE VARIABLES
    void Init(){
               
        U=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        E=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        
    }
    
}

/////////////////////////////////////////////////////////////////////////
//GAUGE TRANSFORMATION AND OPERATIONS TO PERFORM GAUGE TRANSFORMATIONS //
/////////////////////////////////////////////////////////////////////////

namespace GaugeTransformation{
    
    GaugeTransformations *G;

    ///////////////////
    //INITIALIZATION //
    ///////////////////
    
    void Init(){
        
        G=new GaugeTransformations(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);

    }
    
    void SetIdentity(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        #pragma omp parallel for
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    COPY_SUNcMatrix(G->Get(x,y,z),SUNcGroup::UnitMatrix);
                    
                }
            }
        }// END PARALLEL FOR
        
    }
    
    void SetIdentity(){
        
        SetIdentity(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1);
        
    }
    
    void SetRandom(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    RandomNumberGenerator::SUNcMatrix(DOUBLE(2.0)*D_SQRT2*PI,G->Get(x,y,z));
                    
                }
            }
        }
        
    }
    
    void SetRandom(){
        
        SetRandom(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1);
    
    }
    
    void Save(std::string fname){
        Save(fname,G);
    }
    
    //////////////////////////
    //GAUGE TRANSFORMATIONS //
    //////////////////////////
    
    namespace Operations{
        
        //CREATE A COPY OF THE DYNAMICAL FIELDS
        void CopyFields(){
            
            std::memcpy(GaugeFixedVariables::U->Get(0,0,0,0),GLinks::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*GLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            std::memcpy(GaugeFixedVariables::E->Get(0,0,0,0,0),EFields::E->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::VectorSize*EFields::E->Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
            
            
        }
        
        //CREATE A COPY OF THE DYNAMICAL FIELDS
        void SaveFields(){
            
            std::memcpy(GLinks::U->Get(0,0,0,0),GaugeFixedVariables::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*GLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            std::memcpy(EFields::E->Get(0,0,0,0,0),GaugeFixedVariables::E->Get(0,0,0,0,0),Lattice::Dimension*SUNcAlgebra::VectorSize*EFields::E->Volume*sizeof(SU_Nc_ALGEBRA_FORMAT));
            
            
        }
        

        //PERFORM GAUGE TRANSFORMATION OF GAUGE LINKS
        void GaugeTransformLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
            
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        PERFORM_GAUGE_LINK_TRANSFORMATION(x,y,z)
		
                    }
                }
            }// END PARALLEL FOR
    
        }
        
        void GaugeTransformLinks(GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
            
            GaugeTransformLinks(0,UOld->N[0]-1,0,UOld->N[1]-1,0,UOld->N[2]-1,UOld,UNew,G);
            
        }
        
        void GaugeTransformLinks(){
            
            GaugeTransformLinks(GLinks::U,GaugeFixedVariables::U,GaugeTransformation::G);
            
        }
        
        //PERFORM GAUGE TRANSFORMATION OF ELECTRIC FIELDS
        void GaugeTransformElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *EOld,ElectricFields *ENew,GaugeTransformations *G){
            
                #pragma omp parallel for
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            PERFORM_ELECTRIC_FIELD_TRANSFORMATION(x,y,z)
                            
                        }
                    }
                }
            
            
        }
        
        void GaugeTransformElectricFields(ElectricFields *EOld,ElectricFields *ENew,GaugeTransformations *G){
            
            GaugeTransformElectricFields(0,EOld->N[0]-1,0,EOld->N[1]-1,0,EOld->N[2]-1,EOld,ENew,G);
            
        }
        
        void GaugeTransformElectricFields(){
            
            GaugeTransformElectricFields(EFields::E,GaugeFixedVariables::E,GaugeTransformation::G);
            
        }
        
        
    }

    
}

#endif
