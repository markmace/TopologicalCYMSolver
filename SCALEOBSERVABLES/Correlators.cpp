namespace ScaleObservables{
    
    namespace Correlators{
        
        void ComputeCorrelators(INT x1,INT y1,INT z1,INT x2,INT y2,INT z2,GaugeLinks *U,ElectricFields *E,DOUBLE &E1E2,DOUBLE &B1B2){
            
            
            /////////////////////////////
            //GET PARALLEL TRANSPORTER //
            /////////////////////////////
            
            SU_Nc_FUNDAMENTAL_FORMAT U12Transporter[SUNcGroup::MatrixSize];
            
            GetTransporter(x1,y1,z1,x2,y2,z2,U12Transporter,U);
            
            
            ////////////////////////////////
            // COMPUTE LOCAL FIELD VALUES //
            ////////////////////////////////
            
            // SET COMPUTE BUFFERS //
            
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            // SET VALUE ARRAYS //
            DOUBLE ExStart[SUNcAlgebra::VectorSize];
            DOUBLE EyStart[SUNcAlgebra::VectorSize];
            DOUBLE EzStart[SUNcAlgebra::VectorSize];
            
            DOUBLE BxStart[SUNcAlgebra::VectorSize];
            DOUBLE ByStart[SUNcAlgebra::VectorSize];
            DOUBLE BzStart[SUNcAlgebra::VectorSize];
            
            
            // COMPUTE EStart AND BStart //
            COMPUTE_AVG_FIELD_STRENGTH(x1,y1,z1);
            
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                
                ExStart[a]=E0Loc[a]; EyStart[a]=E1Loc[a]; EzStart[a]=E2Loc[a];
                BxStart[a]=B0Loc[a]; ByStart[a]=B1Loc[a]; BzStart[a]=B2Loc[a];

            }
            
            
            // SET VALUE ARRAYS //
            DOUBLE ExEnd[SUNcAlgebra::VectorSize];
            DOUBLE EyEnd[SUNcAlgebra::VectorSize];
            DOUBLE EzEnd[SUNcAlgebra::VectorSize];
            
            DOUBLE BxEnd[SUNcAlgebra::VectorSize];
            DOUBLE ByEnd[SUNcAlgebra::VectorSize];
            DOUBLE BzEnd[SUNcAlgebra::VectorSize];
            
            // COMPUTE EEnd AND BEnd //
            COMPUTE_AVG_FIELD_STRENGTH(x2,y2,z2);
            
            // PARALLEL TRANSPORT TO STARTING POINT  -- COMPUTE U_{x->y}^{ab} B^{b}_{i}(y) //
            
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,E0Loc,ExEnd);
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,E1Loc,EyEnd);
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,E2Loc,EzEnd);

            
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,B0Loc,BxEnd);
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,B1Loc,ByEnd);
            SUNcAlgebra::Operations::AdjointMultiplication(U12Transporter,B2Loc,BzEnd);
            
            // COMPUTE E^{a}_{i}(x) U_{x->y}^{ab} E^{b}_{i}(y) AND B^{a}_{i}(x) U_{x->y}^{ab} B^{b}_{i}(y) //
            
            DOUBLE EXCorr=SUNcAlgebra::Operations::ScalarProduct(ExStart,ExEnd);
            DOUBLE EYCorr=SUNcAlgebra::Operations::ScalarProduct(EyStart,EyEnd);
            DOUBLE EZCorr=SUNcAlgebra::Operations::ScalarProduct(EzStart,EzEnd);
            
            DOUBLE BXCorr=SUNcAlgebra::Operations::ScalarProduct(BxStart,BxEnd);
            DOUBLE BYCorr=SUNcAlgebra::Operations::ScalarProduct(ByStart,ByEnd);
            DOUBLE BZCorr=SUNcAlgebra::Operations::ScalarProduct(BzStart,BzEnd);

            
            E1E2=EXCorr+EYCorr+EZCorr;  B1B2=BXCorr+BYCorr+BZCorr;
            
            
            
        }
        
        void CorrelatorHistogram(std::string fname,INT nPoints,GaugeLinks *U,ElectricFields *E){
            
            // COMMANDLINE READOUT
            std::cerr << "#BEGIN CORRELATOR CALCUALTIONS" << std::endl;
            
            // MAX LENGTH //
            DOUBLE MAX_LENGTH=std::pow((Lattice::N[0]/2*Lattice::a[0])*(Lattice::N[1]/2*Lattice::a[1])*(Lattice::N[2]/2*Lattice::a[2]),1.0/3.0);
            
            INT NumberOfBins=8*INT(MAX_LENGTH);
            
            // CREATE HISTOGRAMS
            Histogram *ECorrelator=new Histogram(0.0,MAX_LENGTH,NumberOfBins);
            Histogram *BCorrelator=new Histogram(0.0,MAX_LENGTH,NumberOfBins);

            for(INT nP=0;nP<nPoints;nP++){
                
                // DETERMINE STARTING POINT //
                INT x1=MOD(INT(RandomNumberGenerator::rng()*(U->N[0])),U->N[0]);
                INT y1=MOD(INT(RandomNumberGenerator::rng()*(U->N[1])),U->N[1]);
                INT z1=MOD(INT(RandomNumberGenerator::rng()*(U->N[2])),U->N[2]);
                
                // DETERMINE END POINT //
                INT x2=MOD(INT(RandomNumberGenerator::rng()*(U->N[0])),U->N[0]);
                INT y2=MOD(INT(RandomNumberGenerator::rng()*(U->N[1])),U->N[1]);
                INT z2=MOD(INT(RandomNumberGenerator::rng()*(U->N[2])),U->N[2]);
                
                // GET DISTANCE //
                DOUBLE L=GetDistance(x1,y1,z1,x2,y2,z2);
                
                // CHECK LENGTH //
                if(L<MAX_LENGTH){
                    
                    // CORRELATORS //
                    DOUBLE B1B2Loc=0; DOUBLE E1E2Loc=0;
                    
                    // COMPUTE CORRELATORS //
                    ComputeCorrelators(x1,y1,z1,x2,y2,z2,U,E,E1E2Loc,B1B2Loc);
                    
                    // COUNT TO HISTOGRAM //
                    ECorrelator->Count(L,E1E2Loc);
                    BCorrelator->Count(L,B1B2Loc);
                    
                }
                
            }
            
            // SET OUTPUT FILE //
            std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,"E",fname,"ID",RandomNumberGenerator::MySEED,".txt");
            std::string BOutputFile=StringManipulation::StringCast(IO::OutputDirectory,"B",fname,"ID",RandomNumberGenerator::MySEED,".txt");
            // SET HEADER FILE //
            std::string HeaderMessage=StringManipulation::StringCast("#L ","CORR ");
            
            // OUTPUT CORRELATORS //
            ECorrelator->Output(HeaderMessage,EOutputFile);
            BCorrelator->Output(HeaderMessage,BOutputFile);

            // CLEAN-UP //
            delete ECorrelator; delete BCorrelator;
            
            // COMMANDLINE READOUT //
            std::cerr << "#FINSIHED CORRELATORS CALCUALTIONS" << std::endl;
        }
        
        // OVERLOAD //
        void CorrelatorHistogram(std::string fname,INT nPoints){
            
            CorrelatorHistogram(fname,nPoints,GLinks::U,EFields::E);
        
        }
        
    }
    
}