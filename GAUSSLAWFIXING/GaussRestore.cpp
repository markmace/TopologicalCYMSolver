namespace GaussLawRestoration{
    
//GAUSS VIOLATION INDEX -- MACRO TO RETURN THE POSITION OF A CERTAIN COMPONENT OF COLORED FIELD
#define GaussViolationIndex(x,y,z,a) ((a)+SUNcAlgebra::VectorSize*(E->Index3D((x),(y),(z))))
    
    // GAUSS VIOLATION ARRAY //
    DOUBLE *GaussViolation;
    
    //GLOBAL MAXIMUM AND AVERAGE VIOLATION OF GAUSS LAW
    DOUBLE GlobalMaxViolation; 
    DOUBLE GlobalAvgViolation;
    
    
    void UpdateDeviation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        
        //MAXIMUM VIOLATION
        DOUBLE MaxViolation=0.0;
        
        //SQR SUM OF ALL VIOLATIONS
        DOUBLE SqrSumViolation=0.0;
        
        
        #pragma omp parallel
        {
            
            //ALLOCATE BUFFERS
            SET_ADJOINT_BUFFERS();
            
            
            // UPDATE ELECTRIC FIELDS //
            #pragma omp for reduction( + : SqrSumViolation) reduction( max : MaxViolation)
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        // GET LOCAL VIOLATION //
                        COMPUTE_NEIGHBORING_ADJOINTS(x,y,z);
                        
                        // COMPUTEVALL COMPONENTS //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            // WRITE TO ARRAY //
                            COMPUTE_GAUSS_VIOLATION(GaussViolation[GaussViolationIndex(x,y,z,a)],x,y,z,a);
                            
                            //UPDATE MAXIMUM
                            MaxViolation=std::max(MaxViolation,GaussViolation[GaussViolationIndex(x,y,z,a)]);
                            
                            //UPDATE SQR SUM
                            SqrSumViolation+=SQR(GaussViolation[GaussViolationIndex(x,y,z,a)]);
                        }
                        
                    }
                }
            }
            
        } // END PARALLEL
        
        //GET GLOBAL RESULTS
        GlobalMaxViolation=MaxViolation; 
        GlobalAvgViolation=sqrt(SqrSumViolation/(SUNcAlgebra::VectorSize*U->Volume));
        
    }
    
    void UpdateDeviation(GaugeLinks *U,ElectricFields *E){

        UpdateDeviation(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);

    }
    
    
    // UPDATE ELECTRIC FIELDS //
    void UpdateElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E,DOUBLE Gamma){
        
        
        #pragma omp parallel
        {
            
            //ALLOCATE BUFFERS
            SU_Nc_ADJOINT_FORMAT U0Adj[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];
            SU_Nc_ADJOINT_FORMAT U1Adj[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];
            SU_Nc_ADJOINT_FORMAT U2Adj[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];
            
            // UPDATE ELECTRIC FIELDS //
            #pragma omp for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        // GET ADJOINTS //
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,0),U0Adj);
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,1),U1Adj);
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,2),U2Adj);
                        
                        // UPDATE ELECTRIC FIELDS //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            for(INT b=0;b<SUNcAlgebra::VectorSize;b++){
                                
                                E->Get(x,y,z,0,a)[0]-=Gamma*(DELTA(a,b)*GaussViolation[GaussViolationIndex(x,y,z,b)]-U0Adj[SUNcGroup::AdjIndex(a,b)]*GaussViolation[GaussViolationIndex(x+1,y,z,b)]);
                                E->Get(x,y,z,1,a)[0]-=Gamma*(DELTA(a,b)*GaussViolation[GaussViolationIndex(x,y,z,b)]-U1Adj[SUNcGroup::AdjIndex(a,b)]*GaussViolation[GaussViolationIndex(x,y+1,z,b)]);
                                E->Get(x,y,z,2,a)[0]-=Gamma*(DELTA(a,b)*GaussViolation[GaussViolationIndex(x,y,z,b)]-U2Adj[SUNcGroup::AdjIndex(a,b)]*GaussViolation[GaussViolationIndex(x,y,z+1,b)]);
                                
                            }
                            
                        }
                        
                        
                    }
                    
                }
            }
        
            
        } // END PARALLEL

        
        
    }
    
    void UpdateElectricFields(GaugeLinks *U,ElectricFields *E,DOUBLE Gamma){
        
        UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E,Gamma);
        
    }

    
    // COMPLETE UPDATE STEP //
    void UpdateStep(GaugeLinks *U,ElectricFields *E,DOUBLE Gamma){
        
        UpdateDeviation(U,E);
        UpdateElectricFields(U,E,Gamma);
        
    }
                                                                       
    
    void Restore(GaugeLinks *U,ElectricFields *E){
        
        // ALLOCATE MEMORY //
        GaussViolation=new DOUBLE[U->Volume*SUNcAlgebra::VectorSize];
        
        // STEP COUNTER //
        INT nSteps=0;
        
        //RESET BREAK CONDITION //
        GlobalMaxViolation=1;
        
        std::cerr << "#GAUSS LAW FIXING" << std::endl;
        
        // ITERATIVE UPDATE //
        while(nSteps<10000 && GlobalMaxViolation>std::pow(10.0,-MAX_DIGITS_PRECISION+1)){
        
            UpdateStep(U,E,7.0/24.0);
            UpdateStep(U,E,7.0/48.0);
            UpdateStep(U,E,7.0/72.0);
            
            nSteps++;
            
            if(nSteps%100==0){
                std::cerr << nSteps << " " << GlobalMaxViolation << " " << GlobalAvgViolation << std::endl;
            }

        
        }
        
        // FREE MEMORY //
        delete[] GaussViolation;
    }
    
    void Restore(){
        Restore(GLinks::U,EFields::E);
    }
    
    
}