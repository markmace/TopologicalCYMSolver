#ifndef __BULKOBSERVABLES__CPP__
#define __BULKOBSERVABLES__CPP__

namespace Observables {
        
    namespace Bulk{
        
        //ELECTRIC AND MAGNETIC FIELD STRENGTH SQUARED
        DOUBLE ESqr[Lattice::Dimension];
        DOUBLE BSqr[Lattice::Dimension];
        
        //CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
        DOUBLE cB[Lattice::Dimension];
        DOUBLE cE[Lattice::Dimension];
        
        //SET CONSTANTS FOR ENERGY MOMENTUM TENSOR
        void SetConstants(GaugeLinks *U){

            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cB[mu]=(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cE[mu]=(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
            
        }
        
        
        // DIAGONAL COMPONENTS OF ENERGY MOMENTUM TENSOR FOR ORIGINAL GaugeLinks
        DOUBLE T00(){
            return 0.5*(cE[0]*ESqr[0]+cE[1]*ESqr[1]+cE[2]*ESqr[2]+cB[0]*BSqr[0]+cB[1]*BSqr[1]+cB[2]*BSqr[2]);
        }
        
        DOUBLE TXX(){
            
            return 0.5*((cE[1]*ESqr[1]+cE[2]*ESqr[2]+cB[1]*BSqr[1]+cB[2]*BSqr[2])-(cE[0]*ESqr[0]+cB[0]*BSqr[0]));
        }
        
        DOUBLE TYY(){
            return 0.5*((cE[0]*ESqr[0]+cE[2]*ESqr[2]+cB[0]*BSqr[0]+cB[2]*BSqr[2])-(cE[1]*ESqr[1]+cB[1]*BSqr[1]));
        }
        
        DOUBLE TZZ(){
            return 0.5*((cE[0]*ESqr[0]+cE[1]*ESqr[1]+cB[0]*BSqr[0]+cB[1]*BSqr[1])-(cE[2]*ESqr[2]+cB[2]*BSqr[2]));
        }
        
        DOUBLE ELECTRIC(){
            return 0.5*(cE[0]*ESqr[0]+cE[1]*ESqr[1]+cE[2]*ESqr[2]);
        }
        
        DOUBLE MAGNETIC(){
            return 0.5*(cB[0]*BSqr[0]+cB[1]*BSqr[1]+cB[2]*BSqr[2]);
        }
        
        
        void UpdateMagnetic(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U){
            
            //SET OBSERVABLES
            DOUBLE B0Sqr=0.0; DOUBLE B1Sqr=0.0; DOUBLE B2Sqr=0.0;
            
            #pragma omp parallel
            {
                //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
                SET_ELEMENTARY_PLAQUETTE_BUFFERS();
                
                //UPDATE ALL GAUGE LINKS
                #pragma omp for reduction( + : B0Sqr,B1Sqr,B2Sqr)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //COMPUTE ELEMENTARY PLAQUETTES AT THIS POINT
                            COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z);
                            
                            //ADD SQUARED FIELD STRENGTH TO LOCAL OBSERVABLES
                            B0Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uyz);
                            B1Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uzx);
                            B2Sqr+=4.0*SUNcGroup::Operations::ReTrIDMinusU(Uxy);
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
            BSqr[0]=B0Sqr/U->Volume; BSqr[1]=B1Sqr/U->Volume; BSqr[2]=B2Sqr/U->Volume;

            
        }
        
        void UpdateMagnetic(GaugeLinks *U){
            
            //SET CONSTANTS TO COMPUTE ENERGY MOMENTUM TENSOR
            SetConstants(U);
            
            UpdateMagnetic(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);

        }
        
        void UpdateElectric(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *E){
            
            //SET OBSERVABLES
            DOUBLE E0Sqr=0.0; DOUBLE E1Sqr=0.0; DOUBLE E2Sqr=0.0;
            
            #pragma omp parallel
            {
                        
                //UPDATE ALL GAUGE LINKS
                #pragma omp for reduction( + : E0Sqr,E1Sqr,E2Sqr)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                                //SUM ALL COLOR COMPONENTS AND ADD TO SQUARED FIELD STRENGTH
                                for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                                    
                                    E0Sqr+=SQR(E->Get(x,y,z,0,a)[0]);
                                    E1Sqr+=SQR(E->Get(x,y,z,1,a)[0]);
                                    E2Sqr+=SQR(E->Get(x,y,z,2,a)[0]);
                                    
                                }                            
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
            ESqr[0]=E0Sqr/E->Volume; ESqr[1]=E1Sqr/E->Volume; ESqr[2]=E2Sqr/E->Volume;
            
        }
        
        
        void Update(GaugeLinks *U,ElectricFields *E){
            
            //SET CONSTANTS TO COMPUTE ENERGY MOMENTUM TENSOR
            SetConstants(U);
            
            //UPDATE INDIVIDUAL CONTRIBUTIONS
            UpdateMagnetic(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
            UpdateElectric(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,E);
            
            
        }
        
        void Update(){
            Update(GLinks::U,EFields::E);
        }
        
        
    }
    
    
}


#endif
