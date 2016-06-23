#ifndef __GAUSSLAW__CPP__
#define __GAUSSLAW__CPP__

namespace Observables {
    
    namespace GaussLaw{
        
        //GLOBAL MAXIMUM AND AVERAGE VIOLATION OF GAUSS LAW
        DOUBLE GlobalMaxViolation;
        DOUBLE GlobalAvgViolation;
        
        //COMPUTES THE VIOLATION OF THE GAUSS LAW CONSTRAINT
        void UpdateViolation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            
            //MAXIMUM VIOLATION
            DOUBLE MaxViolation=0.0;
            
            //SQR SUM OF ALL VIOLATIONS
            DOUBLE SqrSumViolation=0.0;
            
            #pragma omp parallel
            {
                
                //ALLOCATE BUFFERS
                SET_GAUSS_LAW_BUFFERS();
                
                //LOCAL VIOLATION
                DOUBLE LocalViolation[SUNcAlgebra::VectorSize];
                
                //COMPUTE LOCAL VIOLATION AT EACH POINT
                #pragma omp for reduction( + : SqrSumViolation) reduction( max : MaxViolation)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            // GET LOCAL VIOLATION //
                            COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z);
                            
                            //COMPUTE VIOLATION FOR EACH COLOR COMPONENT
                            for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                                
                                //UPDATE MAXIMUM
                                MaxViolation=std::max(MaxViolation,LocalViolation[a]);
                                
                                //UPDATE SQR SUM
                                SqrSumViolation+=SQR(LocalViolation[a]);
                                
                            }
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
            //GET GLOBAL RESULTS
            GlobalMaxViolation=MaxViolation; 
            GlobalAvgViolation=sqrt(SqrSumViolation/(SUNcAlgebra::VectorSize*U->Volume));
            
        }
        
        void CheckViolation(GaugeLinks *U,ElectricFields *E){
            
            UpdateViolation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
            std::cerr << "#GAUSS LAW VIOLATION AT T=" << Dynamics::Time() << " MAX="  << GlobalMaxViolation << " AVG=" << GlobalAvgViolation << std::endl; 
        }
        
        void CheckViolation(){
            CheckViolation(GLinks::U,EFields::E);
        }
        
    }
    
    
}

#endif
