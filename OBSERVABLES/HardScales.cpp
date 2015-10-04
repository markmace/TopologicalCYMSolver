#ifndef __HARD_SCALES__CPP__
#define __HARD_SCALES__CPP__

namespace Observables {
        
    namespace HardScales{
        
        // HARD SCALE TENSOR //
        DOUBLE H[Lattice::Dimension];
        
        // CONSTANTS NEEDED TO COMPUTE HARD SCALES //
        DOUBLE cE[Lattice::Dimension];
        
        //SET CONSTANTS FOR HARD SCALES //
        void SetConstants(GaugeLinks *U){
            
            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cE[mu]=(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
            
        }
        
        
        // CHARACTERISTIC MOMENTUM SCALES //
        DOUBLE LambdaXX(){
            return cE[1]*H[1]+cE[2]*H[2]-cE[0]*H[0];
        }
        
        DOUBLE LambdaYY(){
            return cE[2]*H[2]+cE[0]*H[0]-cE[1]*H[1];
        }
        
        DOUBLE LambdaZZ(){
            return cE[0]*H[0]+cE[1]*H[1]-cE[2]*H[2];
        }

        
        // UPDATE //
        void UpdateHardScales(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U){
            
            // SET OBSERVABLES //
            DOUBLE H0=0.0; DOUBLE H1=0.0; DOUBLE H2=0.0;
            
            #pragma omp parallel
            {
                
                
                // SET DYNAMICAL CONSTANTS -2 \sqrt{-g} a^{3} g^{\mu\alpha} g^{\nu\alpha}/(a_{\mu}^2 a_{\nu}^2) //
                DOUBLE cE[Lattice::Dimension];
                
                DOUBLE gamma=2.0*Dynamics::MetricDeterminant*U->aCube;
                
                for(INT mu=0;mu<Lattice::Dimension;mu++){
                    cE[mu]=Dynamics::gUpMetric[mu]/SQR(U->a[mu]);
                }
                
                //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
                SET_ELEMENTARY_PLAQUETTE_BUFFERS();
                SET_NEIGHBORING_PLAQUETTE_BUFFERS();
                
                //ALLOCATE BUFFERS TO COMPUTE TRACES OF PLAQUETTES
                SET_ELEMENTARY_COLOR_TRACE_BUFFERS();
                SET_NEIGHBORING_COLOR_TRACE_BUFFERS();
                
                //UPDATE ALL GAUGE LINKS
                #pragma omp for reduction( + : H0,H1,H2)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //COMPUTE ELEMENTARY PLAQUETTES AND COLOR TRACES
                            COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z);
                            COMPUTE_ELEMENTARY_COLOR_TRACES();
                            
                            //COMPUTE NEIGHBORING PLAQUETTES AND COLOR TRACES
                            COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z);
                            COMPUTE_NEIGHBORING_COLOR_TRACES();
                            
                            
                            // COMPUTE GAUGE FORCE //
                            for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                                
                                DOUBLE E0Dot=gamma*cE[0]*(cE[1]*(ReTrITaUxy[a]-ReTrITaUxMy[a])-cE[2]*(ReTrITaUzx[a]-ReTrITaUMzx[a]));
                                DOUBLE E1Dot=gamma*cE[1]*(cE[2]*(ReTrITaUyz[a]-ReTrITaUyMz[a])-cE[0]*(ReTrITaUxy[a]-ReTrITaUMxy[a]));
                                DOUBLE E2Dot=gamma*cE[2]*(cE[0]*(ReTrITaUzx[a]-ReTrITaUzMx[a])-cE[1]*(ReTrITaUyz[a]-ReTrITaUMyz[a]));
                                
                                // UPDATE //
                                H0+=SQR(E0Dot); H1+=SQR(E1Dot); H2+=SQR(E2Dot);
                                
                            }
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
            // UPDATE GLOBAL ARRAY //
            H[0]=H0/U->Volume; H[1]=H1/U->Volume; H[2]=H2/U->Volume;
            
        }
        
        void Update(GaugeLinks *U){
            
            //SET CONSTANTS TO COMPUTE ENERGY MOMENTUM TENSOR
            SetConstants(U);
            
            // COMPUTE //
            UpdateHardScales(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);

        }
        
        void Update(){
            Update(GLinks::U);
        }
        
        
    }
    
    
}


#endif
