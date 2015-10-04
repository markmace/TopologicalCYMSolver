#ifndef __JACOBI_GAUGEFIXING__CPP__
#define __JACOBI_GAUGEFIXING__CPP__

namespace JacobiAlgorithm{
    
    DOUBLE alphaStepSize=0.3;
    
    /*MONITORING VARIABLES*/
    DOUBLE MaxStepSize,SqrSumStepSize; 
    DOUBLE MaxDeviation,SqrSumDeviation;  
    DOUBLE MinGaugeFunctional,SumGaugeFunctional;
    
    void UpdateGaugeTransformation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *UGaugeTransformed,GaugeTransformations *G){
        
        //CONVERGENCE MONITORING         
        MaxStepSize=0.0;    MaxDeviation=0.0;   MinGaugeFunctional=(2.0*Nc)*(Dynamics::gUpMetric[0]/SQR(GLinks::U->a[0])+Dynamics::gUpMetric[1]/SQR(GLinks::U->a[1])+Dynamics::gUpMetric[2]/SQR(GLinks::U->a[2])); 
        
        SqrSumDeviation=0.0;    SqrSumStepSize=0.0;     SumGaugeFunctional=0.0;
        
        
        #pragma omp parallel
        {
            //SET BUFFERS
            SET_GAUGE_DEVIATION_BUFFERS();
            SET_GAUGE_UPDATE_BUFFERS();
            
            
            // PERFORM LOCAL UPDATE //
            for(INT z=zLow;z<=zHigh;z++){
                
                #pragma omp for reduction( + : SqrSumDeviation,SqrSumStepSize,SumGaugeFunctional) reduction( max : MaxDeviation,MaxStepSize) reduction( min : MinGaugeFunctional)
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //GET LOCAL GAUGE LINKS
                        GET_LOCAL_GAUGE_VIOLATION(x,y,z);
                        
                        //MONITOR THE QUALITY OF THE GAUGE FIXING 
                        MONITOR_GAUGE_DEVIATION();

                        //COMPUTE JACOBI UPDATE STEP 
                        SUNcAlgebra::Operations::MatrixIExp(-alphaStepSize,DeviationBuffer,wTilde);
                        
                        //MONITOR GAUGE UPDATE
                        MONITOR_GAUGE_UPDATE();
                        
                        //UPDATE GAUGE TRANSFORMATION 
                        SUNcGroup::Operations::DU(wTilde,G->Get(x,y,z),gNew);
                        
                        //WRITE TO GLOBAL ARRAY
                        COPY_SUNcMatrix(G->Get(x,y,z),gNew);
                        
                        
                    }
                    
                }
            }
            
        } // END PARALLEL
        
        
    }
    
    void UpdateGaugeTransformation(GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
        
        //UPDATE EACH CHECKERBOARD PARITY SEPARATELY
        UpdateGaugeTransformation(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,UNew,G);
        GaugeTransformation::Operations::GaugeTransformLinks(UOld,UNew,G);
        
        //SET GLOBAL MONITORING VARIABLES
        GlobalMaxStepSize=MaxStepSize;  GlobalMaxDeviation=MaxDeviation; GlobalMinGaugeFunctional=MinGaugeFunctional;
        
        GlobalAvgStepSize=sqrt(SqrSumStepSize/(GLinks::U->Volume*SUNcAlgebra::VectorSize));   GlobalAvgDeviation=sqrt(SqrSumDeviation/(GLinks::U->Volume*SUNcAlgebra::VectorSize));     GlobalSumGaugeFunctional=SumGaugeFunctional/(GLinks::U->Volume);
        
    }
    
    void UpdateGaugeTransformation(){
        
        UpdateGaugeTransformation(GLinks::U,GaugeFixedVariables::U,GaugeTransformation::G);
        
    }
    
    
}


#endif
