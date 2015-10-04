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
    
    INT NCsOutputFrequency;
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
    
    /////////////////////////////////////////
    //BLOCK GAUGE LINKS AND ELECTRIC FIELDS//
    /////////////////////////////////////////
  
    void RGBlockTypeI(GaugeLinks **UOld,ElectricFields **EOld){
        
        // CREATE NEW OBJECTS //
        GaugeLinks *UNew=new GaugeLinks((*UOld)->N[0]/2,(*UOld)->N[1]/2,(*UOld)->N[2]/2,2*(*UOld)->a[0],2*(*UOld)->a[1],2*(*UOld)->a[2]);
        ElectricFields *ENew=new ElectricFields((*EOld)->N[0]/2,(*EOld)->N[1]/2,(*EOld)->N[2]/2,2*(*EOld)->a[0],2*(*EOld)->a[1],2*(*EOld)->a[2]);
        
        // COMPUTE BLOCKED GAUGE LINKS //
        for(INT z=0;z<=UNew->N[2]-1;z++){
            for(INT y=0;y<=UNew->N[1]-1;y++){
                for(INT x=0;x<=UNew->N[0]-1;x++){
                    
                    //BUFFERS
                    SU_Nc_FUNDAMENTAL_FORMAT UxxP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UyyP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UzzP[SUNcGroup::MatrixSize];
                    
                    SU_Nc_FUNDAMENTAL_FORMAT UzxP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UzxM[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UyxP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UyxM[SUNcGroup::MatrixSize];
                    
                    SU_Nc_FUNDAMENTAL_FORMAT UxyP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UxyM[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UzyP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UzyM[SUNcGroup::MatrixSize];
                    
                    SU_Nc_FUNDAMENTAL_FORMAT UyzP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UyzM[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UxzP[SUNcGroup::MatrixSize];
                    SU_Nc_FUNDAMENTAL_FORMAT UxzM[SUNcGroup::MatrixSize];

                    //SMEARED LINKS FOR BLOCKING
                    SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,0),UxxP);
                    SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,1),UyyP);
                    SUNcGroup::Operations::UU((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,2),UzzP);
                    
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,0),(*UOld)->Get(2*x+1,2*y,2*z+1,0),(*UOld)->Get(2*x+2,2*y,2*z,2),UzxP);
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,1),(*UOld)->Get(2*x+1,2*y+1,2*z,1),(*UOld)->Get(2*x,2*y+2,2*z,0),UxyP);
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,2),(*UOld)->Get(2*x,2*y+1,2*z+1,2),(*UOld)->Get(2*x,2*y,2*z+2,1),UyzP);
                    
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,1),(*UOld)->Get(2*x,2*y+1,2*z,0),(*UOld)->Get(2*x+1,2*y+1,2*z,0),(*UOld)->Get(2*x+2,2*y,2*z,1),UyxP);
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,2),(*UOld)->Get(2*x,2*y,2*z+1,1),(*UOld)->Get(2*x,2*y+1,2*z+1,1),(*UOld)->Get(2*x,2*y+2,2*z,2),UzyP);
                    SUNcGroup::AdvancedOperations::UUUD((*UOld)->Get(2*x,2*y,2*z,0),(*UOld)->Get(2*x+1,2*y,2*z,2),(*UOld)->Get(2*x+1,2*y,2*z+1,2),(*UOld)->Get(2*x,2*y,2*z+2,0),UxzP);

                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y,2*z-1,2),(*UOld)->Get(2*x,2*y,2*z-1,0),(*UOld)->Get(2*x+1,2*y,2*z-1,0),(*UOld)->Get(2*x+2,2*y,2*z-1,2),UzxM);
                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x-1,2*y,2*z,0),(*UOld)->Get(2*x-1,2*y,2*z,1),(*UOld)->Get(2*x-1,2*y+1,2*z,1),(*UOld)->Get(2*x-1,2*y+2,2*z,0),UxyM);
                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y-1,2*z,1),(*UOld)->Get(2*x,2*y-1,2*z,2),(*UOld)->Get(2*x,2*y-1,2*z+1,2),(*UOld)->Get(2*x,2*y-1,2*z+2,1),UyzM);

                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y-1,2*z,1),(*UOld)->Get(2*x,2*y-1,2*z,0),(*UOld)->Get(2*x+1,2*y-1,2*z,0),(*UOld)->Get(2*x+2,2*y-1,2*z,1),UyxM);
                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x,2*y,2*z-1,2),(*UOld)->Get(2*x,2*y,2*z-1,1),(*UOld)->Get(2*x,2*y+1,2*z-1,1),(*UOld)->Get(2*x,2*y+2,2*z-1,2),UzyM);
                    SUNcGroup::AdvancedOperations::DUUU((*UOld)->Get(2*x-1,2*y,2*z,0),(*UOld)->Get(2*x-1,2*y,2*z,2),(*UOld)->Get(2*x-1,2*y,2*z+1,2),(*UOld)->Get(2*x-1,2*y,2*z+2,0),UxzM);
                    
                    //RECONSTRUCT BLOCKED LINKS
                    SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UxxP,UzxP,UzxM,UyxP,UyxM,UNew->Get(x,y,z,0));
                    SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UyyP,UxyP,UxyM,UzyP,UzyM,UNew->Get(x,y,z,1));
                    SUNcGroup::AdvancedOperations::BlockRGTypeIMatrixSum(UzzP,UyzP,UyzM,UxzP,UxzM,UNew->Get(x,y,z,2));
                    
                    // RE-UNITARIZE //
                    SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,0));
                    SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,1));
                    SUNcGroup::Extended::MaxTraceProjection(UNew->Get(x,y,z,2));

                    
                }
            }
        }
        
        // SET ELECTRIC FIELDS TO ZERO //
        ENew->SetZero();
        
        // DELETE OLD OBJECTS //
        delete *UOld; delete *EOld;
        
        // SET POINTER TO BLOCKED FIELDS //
        *UOld=UNew; *EOld=ENew;
        
    }
    
    
    ///////////////////////
    // COOLING PROCEDURE //
    ///////////////////////
    
    void Perform(std::string fname,DOUBLE MaxCoolingTime,INT MaxBlockingNumber,INT Calibrate){
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#COOLING GAUGE LINKS" << std::endl;
        
        //CREATE OUTPUT STREAM FOR LOG FILE //
        std::ofstream OutStream;
        
        std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".log");
        
        OutStream.open(OutputFile.c_str());
        
        //MEASURE BULK OBSERVABLES
        Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
        
        //COMMANDLINE OUTPUT
        DOUBLE VacEst=ChernSimonsNumber::VacuumEstimator::NCS(DynamicFields::U);
        
        OutStream << GradientFlow::CoolingTime() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::Bulk::ELECTRIC() << " " << ChernSimonsNumber::DeltaNCsCooling << " " << VacEst <<  std::endl;
        
        //PERFORM A SEQUENCE OF UPDATES UP TO THE MAXIMUM TIME
        while(GradientFlow::CoolingTime()<MaxCoolingTime){

        
            // PERFORM UPDATE
            if(Calibrate==1){
                GradientFlow::CalibrationUpdate(DynamicFields::U,DynamicFields::E);
            }
            else{
                GradientFlow::StandardUpdate(DynamicFields::U,DynamicFields::E);
            }
            
            if(GradientFlow::CoolingSteps%NCsOutputFrequency==0){
                
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

            if(GradientFlow::CoolingSteps%BlockFrequency==0 && DynamicFields::BlockingNumber<MaxBlockingNumber){
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#BLOCKING LINKS AT COOLING TIME " << GradientFlow::CoolingTime() << std::endl;
                
                //MEASURE BULK OBSERVABLES
                Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
                
                // PERFORM BLOCKING //
                //Block(&DynamicFields::U,&DynamicFields::E);
                //Block(&DynamicFields::U);  Block(&DynamicFields::E);
                RGBlockTypeI(&DynamicFields::U,&DynamicFields::E);

                // INCREASE BLOCKING LEVEL //
                DynamicFields::BlockingNumber++;
                
                //MEASURE BULK OBSERVABLES
                Observables::Bulk::Update(DynamicFields::U,DynamicFields::E);
                
                // COMMANDLINE OUTPUT //
                std::cerr << "#" << DynamicFields::BlockingNumber << " BLOCKING COMPLETED" << std::endl;
                
            }
            
        }
        
        // CLOSE OUTPUT STREAM //
        OutStream.close();
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#COOLING DONE" << std::endl;
        
    }
    
    void CoolNSave(std::string fname,DOUBLE MaxCoolingTime,INT MaxBlockingNumber,GaugeLinks *UHot,GaugeLinks *UCold){
        
        // ALLOCATE //
        Allocate(UHot);
        
        // SET BLOCKING LEVEL TO ZERO //
        DynamicFields::BlockingNumber=0;
        
        // SETUP //
        Setup(UHot);
        
        // RESET GRADIENTFLOW //
        GradientFlow::Reset();
        
        // COOL //
        Perform(fname,MaxCoolingTime,MaxBlockingNumber,0);
        
        // SAVE COOL LINKS //
        SaveCooledLinks(UCold);
        
        // CLEAN-UP //
        Cleanup();
        
    }
    
    void ContinueCooling(std::string fname,DOUBLE MaxCoolingTime,INT MaxBlockingNumber,GaugeLinks *UHot){
        
        // ALLOCATE //
        Allocate(UHot);
        
        // SETUP //
        Setup(UHot);
        
        // COOL //
        Perform(fname,MaxCoolingTime,MaxBlockingNumber,1);
        
        // CLEAN-UP //
        Cleanup();
        
    }
    
}
