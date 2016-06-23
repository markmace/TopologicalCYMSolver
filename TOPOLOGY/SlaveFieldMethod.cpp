#ifndef __SLAVE_FIELD_DYNAMICS__CPP__
#define __SLAVE_FIELD_DYNAMICS__CPP__

namespace ChernSimonsNumber {

    namespace SlaveFieldMethod{
        
        // SLAVE FIELD AT T AND T-DT //
        GaugeTransformations *NewSlaveField;
        GaugeTransformations *OldSlaveField;
        
        // GAUGE TRANSFORMED LINKS AND ELECTRIC FIELDS //
        namespace GaugeTransformedFields{
            
            GaugeLinks *U;
            ElectricFields *E;
            
        }
        
        // GET DYNAMIC LINKS //
        void GetDynamicLinks(){
            
            // PERFORM SLAVE FIELD TRANSFORMATION ON NEW LINKS //
            GaugeTransformation::Operations::GaugeTransformLinks(GLinks::U,GaugeTransformedFields::U,NewSlaveField);
            
        }
        
        ////////////////////
        // INITIALIZATION //
        ////////////////////
        
        void Init(){
            
            NewSlaveField=new GaugeTransformations(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
            OldSlaveField=new GaugeTransformations(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
            
            
            GaugeTransformedFields::U=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
            GaugeTransformedFields::E=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
            
        }
        
        // SET SLAVE FIELD TO IDENTITY //
        void SetIdentity(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeTransformations *G){
            
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        COPY_SUNcMatrix(G->Get(x,y,z),SUNcGroup::UnitMatrix);
                        
                    }
                }
            }// END PARALLEL FOR
            
        }
        
        void SetIdentity(GaugeTransformations *G){
            
            SetIdentity(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,G);
            
        }
        
        // SAVE SLAVE FIELD //
        void Save(std::string fname){
            Save(fname,NewSlaveField);
        }
        
        /////////////////
        // PEAK STRESS //
        /////////////////
        
        DOUBLE PeakStress=0.0;
        
        
        // MEASURE PEAK STRESS //
        void MeasurePeakStress(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,GaugeTransformations *S){
            
            DOUBLE MaxStress=0.0;
            
            #pragma omp parallel
            {
                #pragma omp for reduction(max : MaxStress)
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            // DEFINE LOCAL STRESS //
                            DOUBLE LocalStress=DOUBLE(0.5)*(SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x,y,z,0))+SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x,y,z,1))+SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x,y,z,2))+SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x-1,y,z,0))+SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x,y-1,z,1))+SUNcGroup::Operations::ReTrIDMinusU(GaugeTransformedFields::U->Get(x,y,z-1,2)));
                            
                            // TEST FOR SUPREMUM //
                            MaxStress=std::max(LocalStress,MaxStress);
                            
                        }
                    }
                    
                }
                
            }// END PARALLEL
            
            
            // SET GLOBAL PEAKSTRESS //
            PeakStress=MaxStress;
            
        }
        
        void MeasurePeakStress(GaugeLinks *U,GaugeTransformations *S){
            
            MeasurePeakStress(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,S);
            
        }
        
        //////////////////////////////////////////////////////////////////////////////////////
        // PERFORM GAUGE TRANSFORAMTION DURING SLAVE EVOLUTION PROCESS TO ENSURE SMOOTHNESS //
        //////////////////////////////////////////////////////////////////////////////////////
        
        void PerformTransformation(){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "PERFORMING TRANSFORMATION AT T=" << Dynamics::Time() << " WITH PEAK-STRESS " << PeakStress << std::endl;
            
            // SET FINAL SLAVE TRANSFORMED ELECTRIC FIELDS //
            GaugeTransformation::Operations::GaugeTransformElectricFields(EFields::E,GaugeTransformedFields::E,NewSlaveField);
            
            // COPY ENSLAVED FIELDS TO DYNAMICAL FIELDS //
            Copy(GLinks::U,GaugeTransformedFields::U);
            Copy(EFields::E,GaugeTransformedFields::E);
            
            // SET S(x)=1 //
            SetIdentity(OldSlaveField);
            SetIdentity(NewSlaveField);
            
            // MEASURE PEAK STRESS //
            MeasurePeakStress(GLinks::U,NewSlaveField);
            
            std::cerr << "#PEAK STRESS AT TRANSFORMATION " << PeakStress << std::endl;
            
        }
        
        
        
        namespace QuenchDynamics{
            
            // OUTPUT MONITIOR //
            std::ofstream SlaveFieldMonitor;
            
            /////////////////////////////////////
            // INITIALIZE SLAVE FIELD DYNAMICS //
            /////////////////////////////////////
            
            void Init(std::string fname){
                
                std::cerr << "#INITIALIZING SLAVE FIELD DYNAMICS " << std::endl;
                
                // SET SLAVE FIELD TO IDENTITY //
                SetIdentity(NewSlaveField);
                SetIdentity(OldSlaveField);
                
                // GET DYNAMIC LINKS //
                GetDynamicLinks();
                
                // OPEN OUTPUT FILE //
                SlaveFieldMonitor.open(fname.c_str());
                
                // QUENCH SLAVE FIELD //
                for(INT nSteps=0;nSteps<5000;nSteps++){
                    
                    // UPDATE SLAVE FIELD USING LOS ALAMOS ALGORITHM //
                    CoulombGaugeFixing::LosAlamosAlgorithm::UpdateGaugeTransformation(GLinks::U,GaugeTransformedFields::U,NewSlaveField);
                    
                }
                
                // MEASURE PEAK-STRESS //
                MeasurePeakStress(GLinks::U,NewSlaveField);
                
                // SET COULOMB GAUGE //
                PerformTransformation();
                
            }
            
            /////////////////////////////////////////////////////////
            // SET INITIAL GUESS FOR SLAVE FIELD AT EACH TIME STEP //
            //      S(x+dt)=(S(x,t)S^{D}(x,t-dt))^m S(x,t)         //
            /////////////////////////////////////////////////////////
            
            void SetInitialGuess(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,DOUBLE m,GaugeTransformations *SDOld,GaugeTransformations *SDNew){
                
                #pragma omp parallel for
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            // GET PREVIOUS SLAVE FIELD UPDATE //
                            SU_Nc_FUNDAMENTAL_FORMAT NewSlaveFieldUpdate[SUNcGroup::MatrixSize];
                            
                            SUNcGroup::Operations::UD(SDNew->Get(x,y,z),SDOld->Get(x,y,z),NewSlaveFieldUpdate);
                            
                            // SAVE OLD SLAVE FIELD //
                            COPY_SUNcMatrix(SDOld->Get(x,y,z),SDNew->Get(x,y,z));
                            
                            // CHECK CRITERION FOR MEMORY EFFECTS AND PERFORM UPDATE //
                            DOUBLE MemoryCheck=0.5*SUNcGroup::Operations::ReTrIDMinusU(NewSlaveFieldUpdate);
                            DOUBLE Criteria=SQR(1.0-m);
                            
                            if(MemoryCheck<=Criteria){
                                
                                // COMPUTE UPDATE TO THE POWER m //
                                SUNcGroup::AdvancedOperations::Power(NewSlaveFieldUpdate,NewSlaveFieldUpdate,m);
                                
                                // SET NEW SLAVE FIELD //
                                SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];
                                
                                SUNcGroup::Operations::UU(NewSlaveFieldUpdate,SDNew->Get(x,y,z),Buffer);
                                
                                COPY_SUNcMatrix(SDNew->Get(x,y,z),Buffer);
                                
                            }
                            
                        }
                    }
                }// END PARALLEL
            }
            
            void SetInitialGuess(DOUBLE m,GaugeTransformations *SOld,GaugeTransformations *SNew){
                SetInitialGuess(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,m,SOld,SNew);
            }
            
            ////////////////////////////
            // SLAVE FIELD PARAMETERS //
            ////////////////////////////
            
            // MAXIMAL GAUGE DEVAITION //
            DOUBLE MaxDeviationLimit=1.0;//std::pow(10.0,-2.0);
            
            // STRESS TOLORANCE SMax //
            DOUBLE StressTolerance=1.2;
            
            // BASE FREQUENCY FOR SMOOTHING TRANSFORMATION //
            INT TransformationCounterLimit=5;
            
            // OUTPUT FREQUENCY //
            INT OutputFreq=200;
            
            // COUNTER OF SLAVE FIELD UPDATES //
            INT TransformationCounter=0; 
            
            // OFFSET IN WINDING NUMBER MEASUREMENT //
            INT DeltaNCsOffset=0;
            
            
            // COMPLETE DYNAMICAL UPDATE //
            INT Update(INT NumberOfQuenchSteps){
                
                INT ReturnValue=0;
                
                // CONSTANT FOR REPEATING PREVIOUS STEPS// 
                DOUBLE m=1.0-Dynamics::dTau; 
                
                // EFFECTIVE NUMBER OF STEPS //
                INT EffectiveNumberOfQuenchSteps=NumberOfQuenchSteps;
                
                // ADJUST IF PEAKSTRESS IS LARGE //
                if(PeakStress>StressTolerance){
                    
                    // ADJUST STEP SIZE OF MEMORY UPDATE //
                    m=m*m*m;
                    
                    // ADJUST QUENCHING //
                    EffectiveNumberOfQuenchSteps=3*NumberOfQuenchSteps;
                    
                }
                
                // SET INITIAL GUESS FOR SLAVE FIELD //
                SetInitialGuess(m,OldSlaveField,NewSlaveField);
                
                // GET DYNAMIC LINKS //
                GetDynamicLinks();
                
                INT nSteps=0; CoulombGaugeFixing::LosAlamosAlgorithm::MaxDeviation=1.0;
                
                // QUENCH SLAVE FIELD //
                while(nSteps<EffectiveNumberOfQuenchSteps || CoulombGaugeFixing::LosAlamosAlgorithm::MaxDeviation>MaxDeviationLimit){
                    
                    // UPDATE SLAVE FIELD USING LOS ALAMOS ALGORITHM //
                    CoulombGaugeFixing::LosAlamosAlgorithm::UpdateGaugeTransformation(GLinks::U,GaugeTransformedFields::U,NewSlaveField);
                    
                    // INCREASE STEP COUNTER //
                    nSteps++;
                    
                }
                
                // MEASURE PEAK-STRESS //
                MeasurePeakStress(GLinks::U,NewSlaveField);
                
                // MEASURE WINDING //
                INT DeltaNCs=DeltaNCsOffset+WindingNumber::Measure(NewSlaveField);
                
                // MONITOR SLAVE FIELD //
                if(Dynamics::tSteps%1==0){
                    
                    SlaveFieldMonitor << Dynamics::Time() << " " << DeltaNCs << " " << PeakStress << " " << CoulombGaugeFixing::GlobalMaxDeviation << " " << nSteps << std::endl;
                    
                }
                
                // CHECK PEAK STRESS AND IF SMALL PERFORM TRANSFORMATION //
                if(TransformationCounter>=TransformationCounterLimit && PeakStress<StressTolerance){
                    
                    DeltaNCsOffset=DeltaNCs;
                    
                    TransformationCounter=0;
                    
                    ReturnValue=1;
                    
                }
                
                // INCREASE COUNTER //
                TransformationCounter++;
                
                return ReturnValue;
                
            }
            
        }
        
    }
    
}

#endif
