//LANGEVIN TRANSITION RATE SIMULATION
#define METRIC_FLAG MINKOWSKI_FLAG
#define IC_FLAG THERMAL_FLAG

#include <iostream>
#include <string>

/////////////////////////////////////////
// SU(Nc) GROUP AND ALGEBRA OPERATIONS //
/////////////////////////////////////////

//INCLUDE SU(Nc) ALGEBRA AND GROUP DEFINITION
#include "SUNc/Definitions.cpp"

/////////////////////////////
// RANDOM NUMBER GENERATOR //
/////////////////////////////

//INCLUDE RANDOM NUMBER GENERATOR
#include "MISC/RNG/GSLRandomNumberGenerator.cpp"

///////////////////////////
//  OUTPUT HANDLING      //
///////////////////////////

#include "IO/StringManipulation.cpp"
#include "IO/OutputManagement.cpp"

///////////////////////////
//   INPUT HANDLING      //
///////////////////////////
#if IC_FLAG==LOAD_FLAG
    #include "IO/InputManagement.cpp"
#endif
/////////////////////////////////////////////
// DEFINITION OF LATTICE GRID AND INDEXING //
/////////////////////////////////////////////

#include "LATTICE/3DGrid.cpp"
#include "LATTICE/GaugeLinks.cpp"
#include "LATTICE/ElectricFields.cpp"
#include "LATTICE/GaugeTransformations.cpp"
#include "LATTICE/GLinksEFields.cpp"
#include "LATTICE/Indexing.cpp"
#include "LATTICE/Momenta.cpp"

///////////////////////////
//GENERAL PURPOSE MACROS //
///////////////////////////

//INCLUDE MACROS TO COMPUTE PLAQUETTES

#include "CALC/Plaquettes.cpp"
#include "CALC/AvgFieldStrength.cpp"

//INCLUDE MACROS TO DETERMINE AVERAGE FIELD STRENGTH
#include "CALC/AvgFieldStrength.cpp"

//INCLUDE MACROS TO COMPUTE GAUSS LAW VIOLATION
#include "CALC/GaussViolation.cpp"

//INCLUDE MACROS TO PERFORM GAUGE TRANSFORMATIONS
#include "CALC/GaugeTransformation.cpp"

//INCLUDE UNIMPROVED AND IMPROVED CHERN-SIMONS OPERATORS
#include "CALC/ChernSimonsOperators.cpp"

//INCLUDE ROUTINES TO MEASURE CHERN SIMONS NUMBER DERIVATIVE
#include "TOPOLOGY/ChernSimonsDerivative.cpp"

////////////////////////
// INITIAL CONDITIONS //
////////////////////////

DOUBLE Qs=1.0;
DOUBLE n0=1.0;

///////////////////////
//REAL-TIME DYNAMICS //
///////////////////////

//INCLUDE DYNAMICS
#include "DYNAMICS/Dynamics.cpp"

//////////////////////
//BASIC OBSERVABLES //
//////////////////////

//INCLUDE MEASUREMENT OF BULK OBSERVABLES
#include "OBSERVABLES/BulkObservables.cpp"

//INCLUDE MEASUREMENT OF GAUSS LAW VIOLATION
#include "OBSERVABLES/GaussLawViolation.cpp"

//INCLUDE MEASUREMENT OF UNITARITY VIOLATION
#include "OBSERVABLES/UnitarityViolation.cpp"

//INCLUDE MEASUREMENT OF ENERGY MOMENTUM TENSOR
#include "OBSERVABLES/EnergyMomentumTensor.cpp"

//INCLUDE MEASUREMENT OF HARD SCALES
#include "OBSERVABLES/HardScales.cpp"

//INCLUDE HISTOGRAM
#include "MISC/HISTOGRAM/Histogram.cpp"
#include "MISC/HISTOGRAM/MultiHistogram.cpp"

////////////////////////////////
//GAUSS LAW FIXING PROCEDURES //
////////////////////////////////
#include "GAUSSLAWFIXING/GaussRestore.cpp"

///////////////////////
//INITIAL CONDITIONS //
///////////////////////

#include "INITIALCONDITIONS/SetZero.cpp"
#include "INITIALCONDITIONS/SetRandomMatrices.cpp"
#include "INITIALCONDITIONS/SetQuasiParticles.cpp"

//////////////////////
// THERMAL DYNAMICS //
//////////////////////

#include "DYNAMICS/ThermalDynamics.cpp"

///////////////////////
//REAL-TIME DYNAMICS //
///////////////////////

//INCLUDE DYNAMICS
#include "DYNAMICS/LangevinDynamics.cpp"

////////////////////////////
//GAUGE FIXING PROCEDURES //
////////////////////////////

//INCLUDE GAUGE TRANSFORMATION RULES AND BUFFERS
#include "GAUGETRANSFORMATION/GaugeTransformation.cpp"

//INCLUDE COULOMB GAUGE FIXING PROCEDURE
#include "GAUGETRANSFORMATION/CoulombGaugeFixing.cpp"

//INCLUDE MEASUREMENT OF SPECTRA
#include "OBSERVABLES/Spectra.cpp"

////////////////////////////////////////
//MEASUREMENTS OF CHERN-SIMONS NUMBER //
////////////////////////////////////////

//INCLUDE VACUUM TOPOLOGY ESTIMATOR
#include "TOPOLOGY/VacuumEstimator.cpp"

/////////////////////////////
//YANG-MILLS GRADIENT FLOW //
/////////////////////////////

//GRADIENT FLOW COOLING
#include "GRADIENTFLOW/GradientFlow.cpp"
#include "GRADIENTFLOW/Cooling.cpp"

////////////////////////////////////////
//MEASUREMENTS OF CHERN-SIMONS NUMBER //
////////////////////////////////////////

#include "TOPOLOGY/CoolingMethod.cpp"
#include "TOPOLOGY/WindingNumber.cpp"
#include "TOPOLOGY/SlaveFieldMethod.cpp"
#include "TOPOLOGY/IntegralEstimateMethod.cpp"

//////////////////
//MISCELLANEOUS //
//////////////////

#include "IO/SaveConfiguration.cpp"

#if IC_FLAG==LOAD_FLAG
    #include "IO/LoadConfiguration.cpp"
#endif


//////////
//TESTS //
//////////

#include "MISC/GaugeInvarianceTest.cpp"

///////////////////////////
// SIMULATION PROCEDURE  //
///////////////////////////


namespace Simulation {
    
    // MPI RANDOM NUMBER SEED //
    INT MY_MPI_RNG_SEED;

    /////////////////////////////
    //INITIALIZATION AND EXIT  //
    /////////////////////////////
    
    INT CoolingFrequency,GaugeFixingFrequency,GaugeFixingTotalSteps,ConfigSaveFrequency,CalibFrequency;
    INT NumberOfQuenchSteps;
    
    void Init(){
    
        ///////////////////////////////////
        //SET NCS MEASUREMENT PARAMETERS //
        ///////////////////////////////////
        
        // SET COOLING FREQUENCY //
        CoolingFrequency=20;/*20;*/ ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps=16;//128;//96;
        CalibFrequency=CoolingFrequency*200; ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps=1744;//2336;
        
        
        // SET BLOCKING PARAMETERS //
        Cooling::BlockFrequency=8;//8;//16;
        
        ChernSimonsNumber::CoolingMethod::StandardBlockingLevel=0;
        
        ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel=2;
        
        
        // SLAVE FIELD METHOD //
        NumberOfQuenchSteps=2;
        
        // CONSISTENCY CHECK FOR BLOCKING //
        if((ChernSimonsNumber::CoolingMethod::StandardBlockingLevel*Cooling::BlockFrequency>ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps) || (Cooling::BlockFrequency*ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel>ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps)){
            
            std::cerr << "#ERROR -- INCONSISTENT COOLING PARAMETERS" << std::endl;
            exit(0);
            
        }
        
        // SET GAUGE FIXING AND SPECTRA MEASUREMENT FREQUENCY //
        GaugeFixingFrequency=INT(100.0/(LangevinDynamics::Time())); GaugeFixingTotalSteps=5000;
        
        // SAVE CONFIGURATIONS
        ConfigSaveFrequency=INT(200.0/(LangevinDynamics::Time()));
        
        /////////////////////
        // INITIALIZATIONS //
        /////////////////////
        
        // INITIALIZE 3D LATTICE SETUP //
        Lattice::Init();
        
        // INITIALIZE GAUGE TRANSFORMATIONS //
        GaugeTransformation::Init();
        
        // INITIALIZE GAUGE TRANSFORMED FIELDS //
        GaugeFixedVariables::Init();
        
        // INITIALIZE FOURIER SPACE VARIABLES //
        FourierSpace::Init();
        
        // INITIALIZE COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Init();
        
        // INITIALIZE WINDING NUMBER MEASUREMENT //
        WindingNumber::Init();
        
        // INITIALIZE SLAVE FIELD DYNAMICS
        ChernSimonsNumber::SlaveFieldMethod::Init();
        
        //INITIALIZE LANGEVIN NOISE ELECTRIC FIELD
        LangevinVariables::Init();
        
        
    }
    
    /////////////////////
    //CREATE INFO FILE //
    /////////////////////
    
    void CreateInfoFile(){
        //CREATE INFO FILE //
        std::ofstream OutStream;
        OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"InfoID",MY_MPI_RNG_SEED,".txt").c_str());
        
        //CREATE INFO FILE CONTAINING PARAMETERS
        OutStream << "#INITIAL CONDITIONS" << std::endl;
        OutStream << "#ThermalDynamics::beta= " << ThermalDynamics::beta << std::endl;
        OutStream << "#LangevinDynamics::beta= " << LangevinDynamics::beta << std::endl;
        OutStream << std::endl;
        
        OutStream << "#LATTICE DATA" << std::endl;
        OutStream << "#Nx,Ny,Nz= " << GLinks::U->N[0] << " " << GLinks::U->N[1] << " " << GLinks::U->N[2] << std::endl;
        OutStream << "#ax,ay,az= " << GLinks::U->a[0] << " " << GLinks::U->a[1] << " " << GLinks::U->a[2] << std::endl;
        OutStream << "#dTau= " << Dynamics::dTau << std::endl;
        OutStream << std::endl;
        
        OutStream << "#COOLING PARAMETERS" << std::endl;
        OutStream << "#SqrtDcTauOverSpacing = " << GradientFlow::SqrtDcTauOverSpacing << std::endl;
        OutStream << "#CoolingFrequency = " << CoolingFrequency << std::endl;
        OutStream << "#StandardBlockingLevel= " << ChernSimonsNumber::CoolingMethod::StandardBlockingLevel << std::endl;
        OutStream << "#StandardCoolingMaxSteps= " << ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps << std::endl;
        OutStream << "#BlockFrequency= " << Cooling::BlockFrequency << std::endl;
        OutStream << std::endl;
        
        OutStream << "#CALIBRATION PARAMETERS" << std::endl;
        OutStream << "#CalibFrequency= " << CalibFrequency << std::endl;
        OutStream << "#CalibrationCoolingMaxSteps= " << ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps << std::endl;
        OutStream << "#CalibBlockingLevel= " << ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel << std::endl;
        OutStream << std::endl;
        
        OutStream << "#SLAVE FIELD PARAMETERS" << std::endl;
        OutStream << "#MaxDeviationLimit= " << ChernSimonsNumber::SlaveFieldMethod::QuenchDynamics::MaxDeviationLimit << std::endl;
        OutStream << "#StressTolerance= " << ChernSimonsNumber::SlaveFieldMethod::QuenchDynamics::StressTolerance << std::endl;
        OutStream << "#TransformationCounterLimit= " << ChernSimonsNumber::SlaveFieldMethod::QuenchDynamics::TransformationCounterLimit << std::endl;
        OutStream << "#NumberOfQuenchSteps= " << NumberOfQuenchSteps << std::endl;
        OutStream << std::endl;
        
        OutStream << "#GAUGE FIXING PARAMETERS" << std::endl;
        OutStream << "#GaugeFixingFrequency = " << GaugeFixingFrequency << std::endl;
        OutStream << "#GaugeFixingTotalSteps= " << GaugeFixingTotalSteps << std::endl;
        OutStream << std::endl;
        
        
        // CLOSE OUPUT STREAM //
        OutStream.close();
        

        
        OutStream.close();
        
    }

    ////////////////////////////////////
    // CLASSICAL YANG-MILLS EVOLUTION //
    ////////////////////////////////////
    
    void Evolve(DOUBLE MaxTime){
        
        
        ///////////////////////
        // OUTPUT MANAGEMENT //
        ///////////////////////
        
        std::ofstream EnergyOutStream,CoolNCsOutStream,CalibNCsOutStream;
        
        CoolNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CoolNCsID",MY_MPI_RNG_SEED,".txt").c_str());
        
        CalibNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CalibNCsID",MY_MPI_RNG_SEED,".txt").c_str());
        
        EnergyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"EnergyID",MY_MPI_RNG_SEED,".txt").c_str());

        ///////////////////////////////
        // INITIAL STATE OBSERVABLES //
        ///////////////////////////////
        
        // OPTION TO SAVE INITIAL CONFIGURATION //
        //IO::SaveConfiguration(StringManipulation::StringCast("UOutT",LangevinDynamics::Time()).c_str(),StringManipulation::StringCast("EOutT",LangevinDynamics::Time()).c_str());
        // END OPTION //

        
        // MEASURE BULK OBSERVABLES //
        Observables::Bulk::Update();
        
        // MEASURE HARD SCALES //
        Observables::HardScales::Update();
        
        EnergyOutStream << LangevinDynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
        /*
        // RESET COULOMB GAUGE FIXING ALGORITHMS //
        CoulombGaugeFixing::Reset();
        
        // PERFORM COULOMB GAUGE FIXING //
        CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",LangevinDynamics::Time()).c_str(),GaugeFixingTotalSteps);
        
        // MEASURE SPECTRA
        Observables::Spectra::Compute(StringManipulation::StringCast("SpectraT",LangevinDynamics::Time()).c_str());
         */
        //////////////////////////////////////
        // INITIALIZE TOPOLOGY MEASUREMENTS //
        //////////////////////////////////////
        /*
        // SET SLAVE FIELD INITIALLY TO UNITY AND FIELDS TO GF FIELDS //
        ChernSimonsNumber::SlaveFieldMethod::QuenchDynamics::Init(StringManipulation::StringCast(IO::OutputDirectory,"SlaveID",MY_MPI_RNG_SEED,".txt").c_str());
        
        std::cerr << "#INITIAL PEAKSTRESS AT TIME T=" << Dynamics::Time() << " IS " << ChernSimonsNumber::SlaveFieldMethod::PeakStress <<  std::endl;
        */
        // RESET CHERN SIMONS-MEASUREMENTS //
        ChernSimonsNumber::DeltaNCsRealTime=DOUBLE(0.0);
        ChernSimonsNumber::DeltaNCsCoolRealTime=DOUBLE(0.0);
        ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
        
        // START COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Start();
        
        // CREATE INITIAL COOLING OUTPUT //
        CoolNCsOutStream << LangevinDynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
        /*
        // CALIBRATE //
        ChernSimonsNumber::CoolingMethod::Calibrate();
        
        // CREATE INITIAL CALIBRATION OUTPUT //
        CalibNCsOutStream << LangevinDynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCalibration << " " << ChernSimonsNumber::DeltaNCsPreviousCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimeCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimePreviousCalibration << std::endl;
        */
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        
        while(LangevinDynamics::Time()<MaxTime){

            // UPDATE GAUGE LINKS AND ELECTRIC FIELD VARIABLES //
            LangevinDynamics::Update();
            
            // CHECK BULK OBSERVABLES //
            if(LangevinDynamics::tauLSteps%20==0){
                
                // MEASURE BULK OBSERVABLES //
                Observables::Bulk::Update();
                
                // MEASURE HARD SCALES //
                Observables::HardScales::Update();
                
                EnergyOutStream << LangevinDynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
                
            }
            /*
            // PERFORM GAUGE FIXING //
            if(LangevinDynamics::tauLSteps%GaugeFixingFrequency==0){
                
                // RESET COULOMB GAUGE FIXING ALGORITHMS //
                CoulombGaugeFixing::Reset();
                
                // PERFORM GAUGE FIXING //
                CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",LangevinDynamics::Time()).c_str(),GaugeFixingTotalSteps);
                
                // OPTION TO SAVE GAUGE FIXED FIELDS //
                //GaugeTransformation::Operations::SaveFields();
                // END OPTION //
                
                // MEASURE SPECTRA
                Observables::Spectra::Compute(StringManipulation::StringCast("SpectraT",LangevinDynamics::Time()).c_str());
                
            }
            */
            // UPDATE SLAVE FIELD //
            INT ChangeOfGauge=0;//ChernSimonsNumber::SlaveFieldMethod::QuenchDynamics::Update(NumberOfQuenchSteps);
            
            // COOLED CHERN SIMONS DERIVATIVE //
            if(LangevinDynamics::tauLSteps%CoolingFrequency==0 && ChangeOfGauge==0){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Update(1);
                
                // MEASURE DIFFERNCE IN CHERN SIMONS NUMBER  //
                CoolNCsOutStream << LangevinDynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
            }
            
            // SYNCHRONIZE COOLED IMAGE IN NEW GAUGE //
            if(ChangeOfGauge==1){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Update(1);
                
                // MEASURE DIFFERNCE IN CHERN SIMONS NUMBER  //
                CoolNCsOutStream << LangevinDynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
                
                // CHANGE GAUGE //
                ChernSimonsNumber::SlaveFieldMethod::PerformTransformation();
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Update(0);
                
            }
            
            /*
            // COOLING CALIBRATION //
            if(LangevinDynamics::tauLSteps%CalibFrequency==0){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Calibrate();
                
                // MEASURE OUTPUT //
                CalibNCsOutStream << LangevinDynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCalibration << " " << ChernSimonsNumber::DeltaNCsPreviousCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimeCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimePreviousCalibration << std::endl;
            }
            */
            // OPTION TO SAVE CONFIGURATION //
            /*
             if(LangevinDynamics::tauLSteps%ConfigSaveFrequency==0){
             IO::SaveConfiguration(StringManipulation::StringCast("UOutT",LangevinDynamics::Time()).c_str(),StringManipulation::StringCast("EOutT",LangevinDynamics::Time()).c_str());
             }
             */
            // END OPTION //
            
            // PRECISION CHECKS //
            if(LangevinDynamics::tauLSteps%1000==0){
                
                //CHECK GAUSS LAW VIOLATION
                Observables::GaussLaw::CheckViolation();
                
                //CHECK UNITARITY VIOLATION
                Observables::Unitarity::CheckViolation();
                
            }
            
        }
        
        // CLOSE OUTPUT STREAMS //
        EnergyOutStream.close();   CoolNCsOutStream.close();   CalibNCsOutStream.close();
        
    }
    
    //////////////////////////
    //SIMULATION PROCDEDURE //
    //////////////////////////
    
    void Run(INT MPI_RNG_SEED){
        
        ///////////
        // SETUP //
        ///////////
    
        //SET SEED //
        MY_MPI_RNG_SEED=MPI_RNG_SEED;
        
        //INITIALIZE RANDOM NUMBER GENERATOR //
        RandomNumberGenerator::Init(MY_MPI_RNG_SEED);
    
        //INITIALIZE LANGEVIN DYNAMICS //
        Dynamics::Reset();
        LangevinDynamics::Reset();
        
        //////////////////////
        // CREATE INFO FILE //
        //////////////////////
        
        CreateInfoFile();
        
        //////////////////////////
        // SET INITIAL CONDITON //
        //////////////////////////
        
        #if IC_FLAG==LOAD_FLAG
        IO::LoadConfiguration(GLinks::U,EFields::E);
        #endif
        
        #if IC_FLAG==THERMAL_FLAG
        InitialConditions::SetZero();
        // THERMALIZE
        ThermalDynamics::Thermalize(30,50.0);
        #endif
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        Evolve(4000.0);

    
    }
    
    
}