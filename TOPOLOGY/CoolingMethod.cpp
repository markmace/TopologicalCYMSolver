namespace ChernSimonsNumber{
    
    namespace  CoolingMethod {
        
        INT StandardBlockingLevel; INT CalibrationBlockingLevel;
        
        INT StandardCoolingMaxSteps; INT CalibrationCoolingMaxSteps;
        
        //////////////////////////////////////
        // COOLED GAUGE LINK CONFIGURATIONS //
        //////////////////////////////////////
        
        GaugeLinks *UOld;
        GaugeLinks *UMid;
        GaugeLinks *UNew;
        
        ElectricFields *EMid;
        
        void Init(){
                        
            INT BlockFactor=std::pow(2,StandardBlockingLevel);
            
            UOld=new GaugeLinks(GLinks::U->N[0]/BlockFactor,GLinks::U->N[1]/BlockFactor,GLinks::U->N[2]/BlockFactor,GLinks::U->a[0]*BlockFactor,GLinks::U->a[1]*BlockFactor,GLinks::U->a[2]*BlockFactor);
            UMid=new GaugeLinks(GLinks::U->N[0]/BlockFactor,GLinks::U->N[1]/BlockFactor,GLinks::U->N[2]/BlockFactor,GLinks::U->a[0]*BlockFactor,GLinks::U->a[1]*BlockFactor,GLinks::U->a[2]*BlockFactor);
            UNew=new GaugeLinks(GLinks::U->N[0]/BlockFactor,GLinks::U->N[1]/BlockFactor,GLinks::U->N[2]/BlockFactor,GLinks::U->a[0]*BlockFactor,GLinks::U->a[1]*BlockFactor,GLinks::U->a[2]*BlockFactor);
            
            EMid=new ElectricFields(EFields::E->N[0]/BlockFactor,EFields::E->N[1]/BlockFactor,EFields::E->N[2]/BlockFactor,EFields::E->a[0]*BlockFactor,EFields::E->a[1]*BlockFactor,EFields::E->a[2]*BlockFactor);
            
        }
        
        
        ////////////////////////////
        // PERFORM INTERPOLATION  //
        ////////////////////////////
        
        
        void PerformInterpolation(DOUBLE c,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
            
            DOUBLE cE[Lattice::Dimension];
            
            for(INT mu=0;mu<Lattice::Dimension;mu++){
                
                // NOTE THAT FOR EXPANDING CASE THIS SHOULD REALLY BE THE METRIC IN THE MIDDLE! //
                cE[mu]=(Dynamics::MetricDeterminant*Dynamics::gUpMetric[mu]*UNew->aCube)/(SQR(UNew->a[mu])*Lattice::aScale);
            }
                        
            //UPDATE AT ALL SITES //
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(INT mu=0;mu<Lattice::Dimension;mu++){
                            
                            // GET VALUES //
                            SUNcGroup::AdvancedOperations::GeodesicInterpolation(c,UOld->Get(x,y,z,mu),UNew->Get(x,y,z,mu),UMid->Get(x,y,z,mu),EMid->Get(x,y,z,mu,0));
                            
                            // CONVERT TO STANDARD UNITS //
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                EMid->Get(x,y,z,mu,a)[0]*=cE[mu];
                            }
                            
                        }
                        
                    }
                }
            }
            
        }
        
        void PerformInterpolation(DOUBLE c){
            
            PerformInterpolation(c,0,UNew->N[0]-1,0,UNew->N[1]-1,0,UNew->N[2]-1);
        }
        
        ///////////////////////////////////////
        // PERFORM CHERNS SIMONS MEASUREMENT //
        ///////////////////////////////////////
        
        // IMPROVED SIMPSONS RULE
        DOUBLE GetDeltaNCs(INT NOrder){
            
            // RESET //
            DOUBLE Value=0.0;
            
            // DETERMINE SPACING //
            DOUBLE xSpacing=1.0/DOUBLE(2.0*NOrder);
            
            // GLOBAL WEIGHT //
            DOUBLE GlobalWeight=0.0;
            
            // COMPUTE CONTRIBUTIONS FROM INTERPOLATED POINTS //
            for(INT ip=1;ip<2*NOrder;ip++){
                
                // GET x VALUE //
                DOUBLE xVal=ip*xSpacing;
                
                // GET SIMPSON WEIGHT //
                DOUBLE Weight=2.0+2.0*MOD(ip,2);
                
                // PERFORM INTERPOLATION //
                PerformInterpolation(xVal);
                
                // ADD CONTRIBUTION //
                Value+=Weight*NCsDot(UMid,EMid); GlobalWeight+=Weight;
                
            }
            
            // COMPUTE EMid IF NECESSARY //
            if(NOrder==0){
                PerformInterpolation(0.5);
            }
            
            // COMPUTE LEFT END-POINT CONTRIBUTION //
            Value+=NCsDot(UOld,EMid); GlobalWeight+=1.0;
            
            // COMPUTE RIGHT END-POINT CONTRIBUTION //
            Value+=NCsDot(UNew,EMid);   GlobalWeight+=1.0;
                        
            return  Value/GlobalWeight;
            
        }
        
        //////////////////////////////////
        //  START COOLING METHOD        //
        //////////////////////////////////
        
        void Start(){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#BEGINNING INITIAL COOLING AT T=" << Dynamics::Time() << std::endl;
            
            // PERFORM COOLING //
            Cooling::CoolNSave(StringManipulation::StringCast("CoolingT",Dynamics::Time()).c_str(),StandardCoolingMaxSteps,StandardBlockingLevel,GLinks::U,UNew);
            
            // COPY NEW TO OLD //
            Copy(UOld,UNew);
        }
        
        //////////////////////////////////
        // PERFORM COMPLETE UPDATE STEP //
        //////////////////////////////////
        
        void Update(INT Measure){
            
            // COMMANDLINE OUTPUT //
            //std::cerr << "#COOLING AT T=" << Dynamics::Time() << std::endl;
            
            // RESET CHERN SIMONS NUMBER MONITORING //
            ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
            
            // PERFORM COOLING //
            Cooling::CoolNSave(StringManipulation::StringCast("CoolingT",Dynamics::Time()).c_str(),StandardCoolingMaxSteps,StandardBlockingLevel,GLinks::U,UNew);

            // MEASURE DELTA NCS ALONG COOLING PATH //
            if(Measure==1){
                
                // COMPUTE DIFFERNCE IN CHERN SIMONS NUMBER //
                ChernSimonsNumber::DeltaNCsCoolRealTime+=GetDeltaNCs(1);

            }

            // COPY NEW TO OLD //
            Copy(UOld,UNew);

        }
        
        ///////////////////////////////////////////////
        // CALIBRATE CHERN SIMONS NUMBER MEASUREMENT //
        ///////////////////////////////////////////////
                
        void Calibrate(){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#CALIBRATING COOLING AT T=" << Dynamics::Time() << " Tc=" << GradientFlow::CoolingTime() << std::endl;
        
            
            // RESET CHERN SIMONS NUMBER MONITORING //
            ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
            
            // SAVE PREVIOUS CALIBRATION MEASUREMENTS //
            ChernSimonsNumber::DeltaNCsPreviousCalibration=ChernSimonsNumber::DeltaNCsCalibration;
            ChernSimonsNumber::DeltaNCsRealTimePreviousCalibration=ChernSimonsNumber::DeltaNCsRealTimeCalibration;

            // CALIBRATE BY COOLING TO VACUUM //
            Cooling::ContinueCooling(StringManipulation::StringCast("VacuumCoolingT",Dynamics::Time()).c_str(),StandardCoolingMaxSteps+CalibrationCoolingMaxSteps,CalibrationBlockingLevel,UNew);
            
            // SAVE DELTA NCS ALONG CALIBRATION PATH //
            ChernSimonsNumber::DeltaNCsCalibration=ChernSimonsNumber::DeltaNCsCooling;
            
            // SAVE DELTA NCS ALONG ORIGINAL COOLED REAL-TIME PATH //
            ChernSimonsNumber::DeltaNCsRealTimeCalibration=ChernSimonsNumber::DeltaNCsCoolRealTime;
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#CALIBRATION COMPLETED" << std::endl;

        }
        
    }
    
    
}

