#ifndef __LOS_ALAMOS__CPP__
#define __LOS_ALAMOS__CPP__

namespace LosAlamosAlgorithm{

    DOUBLE OverRelaxationProbability=0.0;
    
    /*MONITORING VARIABLES*/
    DOUBLE MaxStepSize,SqrSumStepSize; 
    DOUBLE MaxDeviation,SqrSumDeviation;  
    DOUBLE MinGaugeFunctional,SumGaugeFunctional;
    
    void UpdateGaugeTransformation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,INT UpdateParity,GaugeLinks *UGaugeTransformed,GaugeTransformations *G){
        
        //CONVERGENCE MONITORING -- RESET ONLY ONCE PER CHECKERBOARD CYCLE
        if(UpdateParity==CHECKERBOARD_EVEN_FLAG){
            
            MaxStepSize=0.0;    MaxDeviation=0.0;   MinGaugeFunctional=(2.0*Nc)*(Dynamics::gUpMetric[0]/SQR(GLinks::U->a[0])+Dynamics::gUpMetric[1]/SQR(GLinks::U->a[1])+Dynamics::gUpMetric[2]/SQR(GLinks::U->a[2])); 
            
            SqrSumDeviation=0.0;    SqrSumStepSize=0.0;     SumGaugeFunctional=0.0;
            
        }
        
        #pragma omp parallel
        {
            
            //SET BUFFERS
            SET_GAUGE_DEVIATION_BUFFERS();
            SET_GAUGE_UPDATE_BUFFERS();
            
            // CHECKERBOARD UPDATE //
            #pragma omp for reduction( + : SqrSumDeviation,SqrSumStepSize,SumGaugeFunctional) reduction( max : MaxDeviation,MaxStepSize) reduction( min : MinGaugeFunctional)
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        if(CHECKERBOARD_PARITY(x,y,z)==UpdateParity || UpdateParity==CHECKERBOARD_ALL_FLAG){
                            
                            //GET LOCAL GAUGE LINKS
                            GET_LOCAL_GAUGE_VIOLATION(x,y,z);
                            
                            //MONITOR THE QUALITY OF THE GAUGE FIXING
                            MONITOR_GAUGE_DEVIATION();
                            
                            //UNITARIZE w(x)
                            SUNcGroup::Extended::MaxTraceProjection(w,wTilde);
                            
                            //MONITOR GAUGE UPDATE
                            MONITOR_GAUGE_UPDATE();
                            
                            //UPDATE GAUGE TRANSFORMATION WITH OVER-RELAXATION
                            if(RandomNumberGenerator::rng()>OverRelaxationProbability){
                                SUNcGroup::Operations::DU(wTilde,G->Get(x,y,z),gNew);
                            }
                            else{
                                SUNcGroup::AdvancedOperations::DDU(wTilde,wTilde,G->Get(x,y,z),gNew);
                            }
                            
                            //WRITE TO GLOBAL ARRAY
                            COPY_SUNcMatrix(G->Get(x,y,z),gNew);
                            
                        }
                    }
                }
            }
            
        } // END PARALLEL
        
        
    }
    
    void UpdateGaugeTransformation(GaugeLinks *UOld,GaugeLinks *UNew,GaugeTransformations *G){
        
        //UPDATE EACH CHECKERBOARD PARITY SEPARATELY
        UpdateGaugeTransformation(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,CHECKERBOARD_EVEN_FLAG,UNew,G);
        GaugeTransformation::Operations::GaugeTransformLinks(UOld,UNew,G);
        
        UpdateGaugeTransformation(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,CHECKERBOARD_ODD_FLAG,UNew,G);
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
