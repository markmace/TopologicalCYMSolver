#ifndef __UNITARITY_VIOLATION__CPP__
#define __UNITARITY_VIOLATION__CPP__

namespace Observables {
    
    namespace Unitarity{
        
        //GLOBAL MAXIMUM AND AVERAGE VIOLATION OF GAUSS LAW
        DOUBLE GlobalMaxViolation;
        DOUBLE GlobalAvgViolation;
        
        //COMPUTES THE VIOLATION OF THE GAUSS LAW CONSTRAINT
        void UpdateViolation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U){
            
            //SET VIOLATIONS BUFFERS
            DOUBLE MaxViolation=0.0; DOUBLE SqrSumViolation=0.0;
            
            #pragma omp parallel
            {
                DOUBLE LocalViolation; 
                
                //COMPUTE LOCAL VIOLATION AT EACH POINT
                #pragma omp for reduction( + : SqrSumViolation) reduction( max : MaxViolation)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //COMPUTE VIOLATION FOR EACH COLOR COMPONENT
                            for(int mu=0;mu<Lattice::Dimension;mu++){
                                
                                //COMPUTE LOCAL VIOLATION
                                LocalViolation=SUNcGroup::Operations::UnitarityNorm(U->Get(x,y,z,mu));
                                
                                //UPDATE MAXIMUM
                                MaxViolation=std::max(MaxViolation,LocalViolation);
                                
                                //UPDATE SQR SUM
                                SqrSumViolation+=SQR(LocalViolation);
                            }
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
            //GET GLOBAL RESULTS
            GlobalMaxViolation=MaxViolation; 
            GlobalAvgViolation=sqrt(SqrSumViolation/(Lattice::Dimension*U->Volume));
            
        }
        
        void CheckViolation(GaugeLinks *U){
            
            UpdateViolation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
            
            std::cerr << "#UNITARITY VIOLATION AT T=" << Dynamics::Time() << " MAX="  << GlobalMaxViolation << " AVG=" << GlobalAvgViolation << std::endl; 
        }
        
        void CheckViolation(DOUBLE CoolTime,GaugeLinks *U){
            
            UpdateViolation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
            
            std::cerr << "#UNITARITY VIOLATION AT T=" << Dynamics::Time() << " Tc=" << CoolTime << " MAX="  << GlobalMaxViolation << " AVG=" << GlobalAvgViolation << std::endl;
        }
        
        void CheckViolation(){
            CheckViolation(GLinks::U);
        }
        
    }
    
    
}

#endif
