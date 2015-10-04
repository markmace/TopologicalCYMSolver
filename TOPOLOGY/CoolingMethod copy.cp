namespace ChernSimonsNumber{
    
    namespace  CoolingMethod {
        
        INT StandardBlockingLevel; INT CalibrationBlockingLevel;
        
        DOUBLE StandardCoolingDepth; DOUBLE CalibrationCoolingDepth;
        
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
        
        
        void PerformInterpolation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
            
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
                            SUNcGroup::AdvancedOperations::GeodesicInterpolation(UOld->Get(x,y,z,mu),UNew->Get(x,y,z,mu),UMid->Get(x,y,z,mu),EMid->Get(x,y,z,mu,0));
                            
                            // CONVERT TO STANDARD UNITS //
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                EMid->Get(x,y,z,mu,a)[0]*=cE[mu];
                            }
                            
                        }
                        
                    }
                }
            }
            
        }
        
        void PerformInterpolation(){
            
            PerformInterpolation(0,UNew->N[0]-1,0,UNew->N[1]-1,0,UNew->N[2]-1);
        }
        
        ///////////////////////////////////////
        // PERFORM CHERNS SIMONS MEASUREMENT //
        ///////////////////////////////////////
        
        DOUBLE GetDeltaNCs(){
            
            return (NCsDot(UOld,EMid)+DOUBLE(4.0)*NCsDot(UMid,EMid)+NCsDot(UNew,EMid))/DOUBLE(6.0); //SIMPSONS RULE
            
        }
        
        
        //////////////////////////////////
        //  START COOLING METHOD        //
        //////////////////////////////////
        
        void Start(){
            
            // PERFORM COOLING //
            Cooling::CoolNSave(StringManipulation::StringCast("CoolingT",Dynamics::Time()).c_str(),StandardCoolingDepth,StandardBlockingLevel,GLinks::U,UNew);
            
            // COPY NEW TO OLD //
            Copy(UOld,UNew);
        }
        
        //////////////////////////////////
        // PERFORM COMPLETE UPDATE STEP //
        //////////////////////////////////
        
        void Update(INT Measure){
            
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#COOLING AT T=" << Dynamics::Time() << std::endl;
            
            // RESET CHERN SIMONS NUMBER MONITORING //
            ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
            
            // PERFORM COOLING //
            Cooling::CoolNSave(StringManipulation::StringCast("CoolingT",Dynamics::Time()).c_str(),StandardCoolingDepth,StandardBlockingLevel,GLinks::U,UNew);

            if(Measure==1)
            // PERFORM INTERPOLATION //
            PerformInterpolation();

            // COMPUTE DIFFERNCE IN CHERN SIMONS NUMBER //
            ChernSimonsNumber::DeltaNCsCoolRealTime+=GetDeltaNCs();

            // COPY NEW TO OLD //
            Copy(UOld,UNew);

        }
        
        ///////////////////////////////////////////////
        // CALIBRATE CHERN SIMONS NUMBER MEASUREMENT //
        ///////////////////////////////////////////////
                
        void Calibrate(){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#CALIBRATING COOLING AT T=" << Dynamics::Time() << std::endl;
            
            // RESET CHERN SIMONS NUMBER MONITORING //
            ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
            
            Cooling::ContinueCooling(StringManipulation::StringCast("VacuumCoolingT",Dynamics::Time()).c_str(),CalibrationCoolingDepth,CalibrationBlockingLevel,UNew);
            
            std::cerr << ChernSimonsNumber::DeltaNCsCooling-ChernSimonsNumber::DeltaNCsPreviousCooling << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << std::endl;
            
            // CALIBRATE NCS MEASUREMENT //
            ChernSimonsNumber::DeltaNCsCoolRealTime=ChernSimonsNumber::DeltaNCsCooling;
            
            
            ChernSimonsNumber::DeltaNCsPreviousCooling=ChernSimonsNumber::DeltaNCsCooling;

            // COMMANDLINE OUTPUT //
            std::cerr << "#CALIBRATION COMPLETED" << std::endl;

        }
        
    }
    
    
}

