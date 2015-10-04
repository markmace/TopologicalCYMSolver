#ifndef __SLAVE_FIELD_DYNAMICS__CPP__
#define __SLAVE_FIELD_DYNAMICS__CPP__

/////////////////
// PEAK STRESS //
/////////////////

DOUBLE PeakStress=0.0;

INT OutputFreq=200;

namespace EnslavedFields{
    
    GaugeLinks *U;
    ElectricFields *E;
    
    //INITIALIZE VARIABLES
    void Init(){
        
        U=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        E=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        
        
    }
    
}

namespace SlaveField{
    
    GaugeTransformations *SDagger;
    GaugeTransformations *SDaggerOld;
    
    //INITIALIZE VARIABLES
    // HERMITIAN CONJUGATE OF SLAVE FIELD
    void Init(){
        
        SDagger=new GaugeTransformations(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        
        SDaggerOld=new GaugeTransformations(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        
        
    }
    // SET SLAVE FIELD TO IDENTITY
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
    
    
    // S^D(x+dt)=(S^D(x,t)S(x,t-dt))^m S^D(x,t)
    void MemoryUpdate(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,DOUBLE m,GaugeTransformations *SDOld,GaugeTransformations *SDNew){
        
        #pragma omp parallel for
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    // GET SLAVE FIELD MOMENTUM //
                    SU_Nc_FUNDAMENTAL_FORMAT SDaggerUpdate[SUNcGroup::MatrixSize];
                    
                    SUNcGroup::Operations::UD(SDNew->Get(x,y,z),SDOld->Get(x,y,z),SDaggerUpdate);
                    
                    // SAVE OLD SLAVE FIELD //
                    COPY_SUNcMatrix(SDOld->Get(x,y,z),SDNew->Get(x,y,z));
                    
                    // CHECK CRITERION FOR MEMORY EFFECTS AND PERFORM UPDATE //
                    DOUBLE MemoryCheck=0.5*SUNcGroup::Operations::ReTrIDMinusU(SDaggerUpdate);
                    DOUBLE Criteria=SQR(1.0-m);
                    
                    if(MemoryCheck<=Criteria){
                        
                        // COMPUTE UPDATE TO THE POWER 1-m //
                        SUNcGroup::AdvancedOperations::Power(SDaggerUpdate,SDaggerUpdate,1.0-m);
                        
                        // SET NEW SLAVE FIELD //
                        SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];
                        
                        SUNcGroup::Operations::UU(SDaggerUpdate,SDNew->Get(x,y,z),Buffer);
                        
                        COPY_SUNcMatrix(SDNew->Get(x,y,z),Buffer);
                        
                    }
                    
                }
            }
        }// END PARALLEL
    }
    
    void MemoryUpdate(DOUBLE m,GaugeTransformations *SOld,GaugeTransformations *SNew){
        MemoryUpdate(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,m,SOld,SNew);
    }

    
    void Save(std::string fname){
        Save(fname,SDagger);
    }
    
    //////////////////////////
    //GAUGE TRANSFORMATIONS //
    //////////////////////////
    
    namespace Operations{
        
        //CREATE A COPY OF THE DYNAMICAL LINKS
        void CopyLinks(){
            
            std::memcpy(EnslavedFields::U->Get(0,0,0,0),GLinks::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*GLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            
        }
        
        //CREATE A COPY OF THE DYNAMICAL LINKS
        void SaveLinks(){
            
            std::memcpy(GLinks::U->Get(0,0,0,0),EnslavedFields::U->Get(0,0,0,0),Lattice::Dimension*SUNcGroup::MatrixSize*GLinks::U->Volume*sizeof(SU_Nc_FUNDAMENTAL_FORMAT));
            
        }
    }
    
}


namespace SlaveFieldDynamics{
    
    ////////////////
    // PARAMETERS //
    ////////////////
    
    DOUBLE MaxDeviationLimit=1.0;//std::pow(10.0,-2.0);
    DOUBLE StressTolerance=1.2;
    INT TransformationCounterLimit=5;
    
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
                        DOUBLE LocalStress=DOUBLE(0.5)*(SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x,y,z,0))+SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x,y,z,1))+SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x,y,z,2))+SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x-1,y,z,0))+SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x,y-1,z,1))+SUNcGroup::Operations::ReTrIDMinusU(EnslavedFields::U->Get(x,y,z-1,2)));
                        
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
    
    // GET DYNAMIC LINKS //
    void GetDynamicLinks(){
        
        // COPY DYNAMIC LINKS TO ENSLAVED LINKS FOR QUENCHING //
        SlaveField::Operations::CopyLinks();
        
        // PERFORM SLAVE FIELD TRANSFORMATION ON NEW LINKS //
        GaugeTransformation::Operations::GaugeTransformLinks(GLinks::U,EnslavedFields::U,SlaveField::SDagger);

    }
    
    // PERFORM FINAL TRANSFORAMTION DURING SLAVE EVOLUTION PROCESS //
    void PerformTransformation(){
       
        // COMMANDLINE OUTPUT //
        std::cerr << "PERFORMING TRANSFORMATION AT T=" << Dynamics::Time() << " WITH PEAK-STRESS " << PeakStress << std::endl;
        
        // SET FINAL SLAVE TRANSFORMED ELECTRIC FIELDS //
        GaugeTransformation::Operations::GaugeTransformElectricFields(EFields::E,EnslavedFields::E,SlaveField::SDagger);
        
        // COPY ENSLAVED FIELDS TO DYNAMICAL FIELDS //
        Copy(GLinks::U,EnslavedFields::U);
        Copy(EFields::E,EnslavedFields::E);
        
        // SET S(x)=1 //
        SlaveField::SetIdentity(SlaveField::SDagger);
        SlaveField::SetIdentity(SlaveField::SDaggerOld);
        
        // MEASURE PEAK STRESS //
        MeasurePeakStress(GLinks::U,SlaveField::SDagger);
        
        std::cerr << "#PEAK STRESS AT TRANSFORMATION " << PeakStress << std::endl;
        
    }
    
    INT TransformationCounter=0; INT DeltaNCsOffset=0;
    
    std::ofstream SlaveFieldMonitor;
    
    void Init(std::string fname){
        
        std::cerr << "#INITIALIZING SLAVE FIELD DYNAMICS " << std::endl;
        
        // SET SLAVE FIELD TO IDENTITY //
        SlaveField::SetIdentity(SlaveField::SDagger);
        SlaveField::SetIdentity(SlaveField::SDaggerOld);
        
        // GET DYNAMIC LINKS //
        GetDynamicLinks();
        
        // OPEN OUTPUT FILE //
        SlaveFieldMonitor.open(fname.c_str());
        
        // QUENCH SLAVE FIELD //
        for(INT nSteps=0;nSteps<5000;nSteps++){
            
            // UPDATE SLAVE FIELD USING LOS ALAMOS ALGORITHM //
            CoulombGaugeFixing::LosAlamosAlgorithm::UpdateGaugeTransformation(GLinks::U,EnslavedFields::U,SlaveField::SDagger);
            
        }
        
        // MEASURE PEAK-STRESS //
        MeasurePeakStress(GLinks::U,SlaveField::SDagger);
        
        // SET COULOMB GAUGE //
        PerformTransformation();
        
    }
    
    INT DissipativeUpdate(INT NumberOfQuenchSteps){
        
        INT ReturnValue=0;
        
        // DISSIPATION CONSTANT AND EFFECTIVE NUMBER OF STEPS //
        DOUBLE m=1.0-Dynamics::dTau; INT EffectiveNumberOfQuenchSteps=NumberOfQuenchSteps;
        
        // ADJUST IF PEAKSTRESS IS LARGE //
        if(PeakStress>StressTolerance){
            
            // ADJUST MEMORY UPDATE //
            m=m*m*m;
            
            // ADJUST QUENCHING //
            EffectiveNumberOfQuenchSteps=3*NumberOfQuenchSteps;
            
        }
        
        // PERFORM MEMORY UPDATE //
        SlaveField::MemoryUpdate(m,SlaveField::SDaggerOld,SlaveField::SDagger);
        
        // GET DYNAMIC LINKS //
        GetDynamicLinks();
        
        INT nSteps=0; CoulombGaugeFixing::LosAlamosAlgorithm::MaxDeviation=1.0;
        
        // QUENCH SLAVE FIELD //
        while(nSteps<EffectiveNumberOfQuenchSteps || CoulombGaugeFixing::LosAlamosAlgorithm::MaxDeviation>MaxDeviationLimit){
            
            // UPDATE SLAVE FIELD USING LOS ALAMOS ALGORITHM //
            CoulombGaugeFixing::LosAlamosAlgorithm::UpdateGaugeTransformation(GLinks::U,EnslavedFields::U,SlaveField::SDagger);
            
            // INCREASE STEP COUNTER //
            nSteps++;
            
        }
        
        // MEASURE PEAK-STRESS //
        MeasurePeakStress(GLinks::U,SlaveField::SDagger);
        
        // MEASURE WINDING //
        INT DeltaNCs=DeltaNCsOffset+ChernSimonsNumber::SlaveFieldMethod::NCS(SlaveField::SDagger);
        
        // MONITOR SLAVE FIELD //
        if(Dynamics::tSteps%1==0){
            
            SlaveFieldMonitor << Dynamics::Time() <<  " " << DeltaNCs << " " << PeakStress << " " << CoulombGaugeFixing::GlobalMaxDeviation << " " << nSteps << std::endl;
            
        }
        
        // CHECK PEAK STRESS AND IF SMALL PERFORM TRANSFORMATION //
        if(TransformationCounter>=TransformationCounterLimit && PeakStress<StressTolerance){
            
            //PerformTransformation();
            
            DeltaNCsOffset=DeltaNCs;
            
            TransformationCounter=0;
            
            ReturnValue=1;
            
        }
        
        // INCREASE COUNTER //
        TransformationCounter++;
        
        return ReturnValue;
        
    }
    
}

#endif
