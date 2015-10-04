#ifndef __TIME_WILSON_LOOP__CPP__
#define __TIME_WILSON_LOOP__CPP__

namespace Observables{
    
    namespace TimeWilsonLoop{
        
        ////////////////
        // PARAMETERS //
        ////////////////
        
        INT NumberOfSamples=40000;
        
        INT tStepRef;
        
        INT TimeRef;
    
        //////////////////////////
        //   REFERENCE FIELDS   //
        //////////////////////////
        namespace ReferenceFields{
            
            GaugeLinks *BlockedU;

        }
        /////////////////////////////////////////////
        // LINKS FOR WILSON LOOPS -- TO BE BLOCKED //
        /////////////////////////////////////////////
        GaugeLinks *BlockedU;
        
        //INITIALIZE VARIABLE
        void Init(){
            
            ReferenceFields::BlockedU=new GaugeLinks(GLinks::U->N[0]/2,GLinks::U->N[1]/2,GLinks::U->N[2]/2,2*GLinks::U->a[0],2*GLinks::U->a[1],2*GLinks::U->a[2]);
            
            BlockedU=new GaugeLinks(GLinks::U->N[0]/2,GLinks::U->N[1]/2,GLinks::U->N[2]/2,2*GLinks::U->a[0],2*GLinks::U->a[1],2*GLinks::U->a[2]);
            
        }

        // BEGIN MEASUREMENT OF SPACE/TIME WILSON LOOP BY SETTING REFERENCE LINKS
        void Begin(GaugeLinks *U){
            
            //SET CONSTANTS AND COPY GAUGE LINKS FROM DYNAMICS
            TimeRef=Dynamics::Time();
            tStepRef=Dynamics::tSteps;
            
            std::cerr << "#BEGIN BLOCKING REFERENCE FIELDS" << std::endl;

            // BLOCK REFERENCE FIELDS //
            CopyWithBlock(&U,ReferenceFields::BlockedU);
            
            std::cerr << "#DONE BLOCKING REFERENCE FIELDS" << std::endl;
        
        }
        
        void Begin(){
        
            Begin(GLinks::U);

        }

        // TRAVERSE THROUGH LATTICE
        void GoXPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UU(Ux,U->Get(x+Var,y,z,0),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UU(UxTmp,U->Get(x+Var,y,z,0),Ux);
                }
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
                
            }
            
        }
        
        void GoXMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UD(Ux,U->Get(x-Var-1,y,z,0),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UD(UxTmp,U->Get(x-Var-1,y,z,0),Ux);
                }
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
            }
            
        }
        
        void GoYPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UU(Ux,U->Get(x,y+Var,z,1),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UU(UxTmp,U->Get(x,y+Var,z,1),Ux);
                }
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
            }
            
        }
        
        void GoYMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UD(Ux,U->Get(x,y-Var-1,z,1),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UD(UxTmp,U->Get(x,y-Var-1,z,1),Ux);
                }
                
                
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
            }
            
        }
        
        void GoZPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UU(Ux,U->Get(x,y,z+Var,2),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UU(UxTmp,U->Get(x,y,z+Var,2),Ux);
                }
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
            }
            
        }
        
        void GoZMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
            
            SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
            
            for(INT Var=0;Var<L;Var++){
                
                if(Var%2==0){
                    SUNcGroup::Operations::UD(Ux,U->Get(x,y,z-Var-1,2),UxTmp);
                }
                else{
                    SUNcGroup::Operations::UD(UxTmp,U->Get(x,y,z-Var-1,2),Ux);
                }
            }
            
            if((L+1)%2==0){
                COPY_SUNcMatrix(Ux,UxTmp);
            }
            
        }
        
        // MEASURES CURRENT CONFIGURATION IN COMPARISON TO REFERENCE TIMES CONFIGURATION
        void Measure(INT L,INT NSamples,GaugeLinks *URef,GaugeLinks *U,DOUBLE &UxTAvg,DOUBLE &UyTAvg,DOUBLE &UzTAvg){
            
            // RESET AVERAGE //
            UxTAvg=DOUBLE(0.0); UyTAvg=DOUBLE(0.0); UzTAvg=DOUBLE(0.0);
            
            #pragma omp parallel
            {
                
                // BUFFERS FOR WILSON LOOP IN SPACE AND TIME
                SU_Nc_FUNDAMENTAL_FORMAT UxT[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UyT[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT UzT[SUNcGroup::MatrixSize];
                
                #pragma omp for reduction( + : UxTAvg,UyTAvg,UzTAvg)
                for(INT n=0;n<NSamples;n++){
                    
                    // DETERMINE PARTITION //
                    INT ID=omp_get_thread_num(); INT NumberOfThreads=omp_get_num_threads();
                    
                    INT SizePerProcess=(U->N[0]*U->N[1]*U->N[2])/NumberOfThreads;  INT Remainder=(U->N[0]*U->N[1]*U->N[2])-NumberOfThreads*SizePerProcess;
                    
                    INT MaxIndex=(ID+1)*(U->N[0]*U->N[1]*U->N[2])/NumberOfThreads;
                    
                    if(ID==0){
                        MaxIndex+=Remainder;
                    }
                    
                    // SET RANDOM BASE POINT //
                    INT Index=MOD(INT(RandomNumberGenerator::rng()*MaxIndex),(U->N[0]*U->N[1]*U->N[2]));
                    
                    INT x,y,z;  U->GetPosition(Index,x,y,z);
                    
                    // CHECK CONSISTENCY //
                    if(U->Index3D(x,y,z)!=Index){
                        
                        std::cerr << "## ERROR -- INDEX IS WRONG" << std::endl;
                        
                        exit(0);
                    }
                    /*
                    // SET RANDOM  BASE POINT //
                    INT x=MOD(INT(RandomNumberGenerator::rng()*(U->N[0])),U->N[0]);
                    INT y=MOD(INT(RandomNumberGenerator::rng()*(U->N[1])),U->N[1]);
                    INT z=MOD(INT(RandomNumberGenerator::rng()*(U->N[2])),U->N[2]);
                    */
                    
                    // SET LINKS TO UNITY //
                    COPY_SUNcMatrix(UxT,SUNcGroup::UnitMatrix);
                    COPY_SUNcMatrix(UyT,SUNcGroup::UnitMatrix);
                    COPY_SUNcMatrix(UzT,SUNcGroup::UnitMatrix);
                    
                    // COMPUTE SEGMENTS //
                    GoXPlus(x,y,z,L,UxT,URef);
                    GoXMinus(x+L,y,z,L,UxT,U);
                    
                    GoYPlus(x,y,z,L,UyT,URef);
                    GoYMinus(x,y+L,z,L,UyT,U);
                    
                    GoZPlus(x,y,z,L,UzT,URef);
                    GoZMinus(x,y,z+L,L,UzT,U);
                    
                    // COMPUTE TRACE //
                    UxTAvg+=SUNcGroup::Operations::ReTr(UxT);
                    UyTAvg+=SUNcGroup::Operations::ReTr(UyT);
                    UzTAvg+=SUNcGroup::Operations::ReTr(UzT);
                    
                }
            } // END PARALLEL
            
            // NORMALIZE //
            UxTAvg=UxTAvg/DOUBLE(Nc*NSamples);
            UyTAvg=UyTAvg/DOUBLE(Nc*NSamples);
            UzTAvg=UzTAvg/DOUBLE(Nc*NSamples);
            
        }

        void WilsonLoopHistogram(std::ofstream *OutStream,GaugeLinks *U){
            
            // COMMANDLINE READOUT
            std::cerr << "#BEGIN TIME WILSON LOOP CALCUALTIONS" << std::endl;
            
            // BLOCK LINKS //
            CopyWithBlock(&U,BlockedU);
    
            // GET DELTA T //
            DOUBLE dT=Dynamics::Time()-TimeRef;
                        
            // STREAM OUT FOR ALL SPATIAL LENGTHS POSSIBLE
            for (INT L=1;L<(BlockedU->N[0])/2;L++){
                
                DOUBLE UxTAvg,UyTAvg,UzTAvg;
                
                Measure(L,NumberOfSamples,ReferenceFields::BlockedU,BlockedU,UxTAvg,UyTAvg,UzTAvg);
                    
                *OutStream << dT << " " << Qs*L*BlockedU->a[0] << " " << UxTAvg << " " << UyTAvg << " " << UzTAvg << std::endl;
                
            }
            
            // NEWLINE OUTPUT STREAM //
            *OutStream << std::endl;
            
            // COMMANDLINE READOUT
            std::cerr << "#FINSIHED TIME WILSON LOOP CALCUALTIONS" << std::endl;
            
        }
        
        void WilsonLoopHistogram(std::ofstream *OutStream){
            
            WilsonLoopHistogram(OutStream,GLinks::U);
            
        }
        
    }
    
    
}

#endif
