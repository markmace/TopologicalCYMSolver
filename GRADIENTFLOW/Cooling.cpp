namespace Cooling{
    
    
    //////////////////////////////
    //   GRADIENT FLOW FIELDS   //
    //////////////////////////////
    
    namespace DynamicFields{
        
        GaugeLinks *U;
        
        ElectricFields *E;
        
        INT BlockingNumber;
        
    }
    ////////////////
    // PARAMETERS //
    ////////////////
    
    INT NCsOutputFrequency=24;
    INT BlockFrequency;
    
    
    ///////////
    // SETUP //
    ///////////
    
    void Allocate(GaugeLinks *UHot){
        
        // ALLOCATE //
        DynamicFields::U=new GaugeLinks(UHot->N[0],UHot->N[1],UHot->N[2],UHot->a[0],UHot->a[1],UHot->a[2]);
        DynamicFields::E=new ElectricFields(UHot->N[0],UHot->N[1],UHot->N[2],UHot->a[0],UHot->a[1],UHot->a[2]);
        
    }
    
    void Setup(GaugeLinks *UHot){
        
        // COPY HOT FIELDS //
        Copy(DynamicFields::U,UHot);
        
        // SET ELECTRIC FIELDS TO ZERO //
        DynamicFields::E->SetZero();
        
    }
    
    void Cleanup(){
        
        delete DynamicFields::U;
        delete DynamicFields::E;
        
    }
    
    
    ////////////////////////
    // SAVE COOLED FIELDS //
    ////////////////////////
    
    void SaveCooledLinks(GaugeLinks *UCold){
        
        Copy(UCold,DynamicFields::U);
        
    }
    
    ///////////////////////
    // COOLING PROCEDURE //
    ///////////////////////
    
    void Perform(std::string fname,INT MaxCoolingSteps,INT MaxBlockingNumber,INT Calibrate){
        
        // COMMANDLINE OUTPUT //
        //std::cerr << "#COOLING GAUGE LINKS" << std::endl;
        
        //CREATE OUTPUT STREAM FOR LOG FILE //
        std::ofstream OutStream;
        
        if(Calibrate==1){
            
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".log");
            
            OutStream.open(OutputFile.c_str());
        }
        
        //MEASURE BULK OBSERVABLES
        Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
        
        //COMMANDLINE OUTPUT
        DOUBLE VacEst=ChernSimonsNumber::VacuumEstimator::NCS(DynamicFields::U);
        
        OutStream << GradientFlow::CoolingTime() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::Bulk::ELECTRIC() << " " << ChernSimonsNumber::DeltaNCsCooling << " " << VacEst <<  std::endl;
        
        //PERFORM A SEQUENCE OF UPDATES UP TO THE MAXIMUM TIME
        while(GradientFlow::CoolingSteps<MaxCoolingSteps){
            
            // PERFORM UPDATE
            if(Calibrate==1){
                GradientFlow::CalibrationUpdate(DynamicFields::U,DynamicFields::E);
            }
            else{
                GradientFlow::StandardUpdate(DynamicFields::U,DynamicFields::E);
            }
            // CREATE OUTPUT //
            if(Calibrate==1){
                
                //MEASURE BULK OBSERVABLES
                Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
                
                //COMMANDLINE OUTPUT
                DOUBLE VacEst=ChernSimonsNumber::VacuumEstimator::NCS(DynamicFields::U);
                
                OutStream << GradientFlow::CoolingTime() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::Bulk::ELECTRIC() << " " << ChernSimonsNumber::DeltaNCsCooling << " " << VacEst <<  std::endl;
                
            }
            
            if(GradientFlow::CoolingSteps%100==0){
                
                DOUBLE CoolTime=GradientFlow::CoolingTime();
                
                //CHECK UNITARITY VIOLATION
                Observables::Unitarity::CheckViolation(CoolTime,DynamicFields::U);
            }
            
            if(GradientFlow::CoolingSteps%(BlockFrequency*INT(std::pow(8,DynamicFields::BlockingNumber)))==0 && DynamicFields::BlockingNumber<MaxBlockingNumber){
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#BLOCKING LINKS AT COOLING TIME " << GradientFlow::CoolingTime() << std::endl;
                
                //MEASURE BULK OBSERVABLES
                Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
                
                // PERFORM BLOCKING //
                Block(&DynamicFields::U);  Block(&DynamicFields::E);
                
                // INCREASE BLOCKING LEVEL //
                DynamicFields::BlockingNumber++;
                
                //MEASURE BULK OBSERVABLES
                Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#" << DynamicFields::BlockingNumber << " BLOCKING COMPLETED" << std::endl;
                
            }
            
        }
        // CLOSE OUTPUT STREAM //
        if(Calibrate==1){
            OutStream.close();
        }
        
        // COMMANDLINE OUTPUT //
        //std::cerr << "#COOLING DONE" << std::endl;
        
    }
    
    void CoolNSave(std::string fname,DOUBLE MaxCoolingSteps,INT MaxBlockingNumber,GaugeLinks *UHot,GaugeLinks *UCold){
        
        // ALLOCATE //
        Allocate(UHot);
        
        // SET BLOCKING LEVEL TO ZERO //
        DynamicFields::BlockingNumber=0;
        
        // SETUP //
        Setup(UHot);
        
        // RESET GRADIENTFLOW //
        GradientFlow::Reset();
        
        // COOL //
        Perform(fname,MaxCoolingSteps,MaxBlockingNumber,0);
        
        // SAVE COOL LINKS //
        SaveCooledLinks(UCold);
    
        // CLEAN-UP //
        Cleanup();
        
    }
    
    void ContinueCooling(std::string fname,DOUBLE MaxCoolingSteps,INT MaxBlockingNumber,GaugeLinks *UHot){
        
        // ALLOCATE //
        Allocate(UHot);
        
        // SETUP //
        Setup(UHot);
        
        // COOL //
        Perform(fname,MaxCoolingSteps,MaxBlockingNumber,1);
        
        // CLEAN-UP //
        Cleanup();
        
    }
    
}
