#ifndef __SPATIAL_WILSON_LOOP__CPP__
#define __SPATIAL_WILSON_LOOP__CPP__

namespace ScaleObservables {

    
    namespace SpatialWilsonLoop{
        
        /////////////////////////////////////////////
        // LINKS FOR WILSON LOOPS -- TO BE BLOCKED //
        /////////////////////////////////////////////
        GaugeLinks *BlockedU;
        
        void Init(){
            
            // ALLOCATE //
            BlockedU=new GaugeLinks(GLinks::U->N[0]/2,GLinks::U->N[1]/2,GLinks::U->N[2]/2,2*GLinks::U->a[0],2*GLinks::U->a[1],2*GLinks::U->a[2]);
            
        }
        
        void Cleanup(){
            
            delete BlockedU;
            
        }
    
        // COMPUTES Tr(WILSON LOOP) OF A GIVEN SIZE AT A GIVEN POINT
        void CalculateWilsonLoop(INT L,INT NSamples,GaugeLinks *U,DOUBLE &UxyAvg,DOUBLE &UyzAvg,DOUBLE &UzxAvg){
            
            // RESET AVERAGE //
            UxyAvg=DOUBLE(0.0); UyzAvg=DOUBLE(0.0); UzxAvg=DOUBLE(0.0);
            #pragma omp parallel
            {
                
                // BUFFERS FOR WILSON LOOP SIDES IN X-Y PLANE
                SU_Nc_FUNDAMENTAL_FORMAT Uxy[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT Uyz[SUNcGroup::MatrixSize];
                SU_Nc_FUNDAMENTAL_FORMAT Uzx[SUNcGroup::MatrixSize];
                
                
                #pragma omp for reduction( + : UxyAvg,UyzAvg,UzxAvg)
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
                    COPY_SUNcMatrix(Uxy,SUNcGroup::UnitMatrix);
                    COPY_SUNcMatrix(Uyz,SUNcGroup::UnitMatrix);
                    COPY_SUNcMatrix(Uzx,SUNcGroup::UnitMatrix);
                    
                    // COMPUTE SEGMENTS //
                    GoXPlus(x,y,z,L,Uxy,U);
                    GoYPlus(x+L,y,z,L,Uxy,U);
                    GoXMinus(x+L,y+L,z,L,Uxy,U);
                    GoYMinus(x,y+L,z,L,Uxy,U);
                    
                    GoYPlus(x,y,z,L,Uyz,U);
                    GoZPlus(x,y+L,z,L,Uyz,U);
                    GoYMinus(x,y+L,z+L,L,Uyz,U);
                    GoZMinus(x,y,z+L,L,Uyz,U);
                    
                    GoZPlus(x,y,z,L,Uzx,U);
                    GoXPlus(x,y,z+L,L,Uzx,U);
                    GoZMinus(x+L,y,z+L,L,Uzx,U);
                    GoXMinus(x+L,y,z,L,Uzx,U);
                    
                    // COMPUTE TRACE //
                    UxyAvg+=SUNcGroup::Operations::ReTr(Uxy);
                    UyzAvg+=SUNcGroup::Operations::ReTr(Uyz);
                    UzxAvg+=SUNcGroup::Operations::ReTr(Uzx);
                    
                }
            } // END PARALLEL
            
            // NORMALIZE //
            UxyAvg=UxyAvg/DOUBLE(Nc*NSamples);
            UyzAvg=UyzAvg/DOUBLE(Nc*NSamples);
            UzxAvg=UzxAvg/DOUBLE(Nc*NSamples);
            
        }
        
          // BLOCKED //
        void BlockedWilsonLoopHistogram(std::string fname,GaugeLinks *U){
            
            // COMMANDLINE READOUT
            std::cerr << "#BEGIN WILSON LOOP CALCUALTIONS" << std::endl;
            
            // COPY DYNAMIC LINKS TO BLOCKED LINKS //
            // OPTION TO BLOCK //
            CopyWithBlock(&U,BlockedU);
            // END OPTION //
            
            // SETUP OUTPUT //
            std::ofstream OutStream;
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            OutStream.open(OutputFile.c_str());
            
            // STREAM OUT FOR ALL SPATIAL LENGTHS POSSIBLE //
            for (INT L=1;L<(BlockedU->N[0])/2;L++){
                
                DOUBLE UxyAvg,UyzAvg,UzxAvg;
                
                CalculateWilsonLoop(L,40000,BlockedU,UxyAvg,UyzAvg,UzxAvg);
                
                OutStream << L*BlockedU->a[0]*Qs << " " << UxyAvg  << " " << UyzAvg  << " " << UzxAvg << std::endl;
                
            }
            
            // CLOSE OUTPUT STREAM //
            OutStream.close();
            
            
            
            // COMMANDLINE READOUT //
            std::cerr << "#FINSIHED WILSON LOOP CALCUALTIONS" << std::endl;
        }
        void BlockedWilsonLoopHistogram(std::string fname){
            
            BlockedWilsonLoopHistogram(fname,GLinks::U);
            
        }
        
        // UNBLOCKED //
        void WilsonLoopHistogram(std::string fname,GaugeLinks *U){
            
            // COMMANDLINE READOUT
            std::cerr << "#BEGIN WILSON LOOP CALCUALTIONS" << std::endl;

            // SETUP OUTPUT //
            std::ofstream OutStream;
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            OutStream.open(OutputFile.c_str());
            
            // STREAM OUT FOR ALL SPATIAL LENGTHS POSSIBLE //
            for (INT L=1;L<(U->N[0])/2;L++){
                
                DOUBLE UxyAvg,UyzAvg,UzxAvg;
                
                CalculateWilsonLoop(L,40000,U,UxyAvg,UyzAvg,UzxAvg);
                
                OutStream << L*U->a[0]*Qs << " " << UxyAvg  << " " << UyzAvg  << " " << UzxAvg << std::endl;
                
            }
            
            // CLOSE OUTPUT STREAM //
            OutStream.close();
            
            // COMMANDLINE READOUT //
            std::cerr << "#FINSIHED WILSON LOOP CALCUALTIONS" << std::endl;
        }
        void WilsonLoopHistogram(std::string fname){
            
            WilsonLoopHistogram(fname,GLinks::U);
            
        }
        
    }
    
}

#endif
