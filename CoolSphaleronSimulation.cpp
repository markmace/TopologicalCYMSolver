//SPHALERON TRANSITION RATE SIMULATION
#define METRIC_FLAG MINKOWSKI_FLAG
#define INPUT_FLAG NO_FLAG
#define NAIVENCS_FLAG NO_FLAG
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
#if INPUT_FLAG==YES_FLAG
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

//INCLUDE SOERENS IMPROVED OPERATORS
#include "CALC/ImprovedOperators.cpp"

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
#include "INITIALCONDITIONS/RestartQuasiParticles.cpp"

//////////////////////
// THERMAL DYNAMICS //
//////////////////////

//INCLUDE DYNAMICS
#if IC_FLAG==THERMAL_FLAG
#include "DYNAMICS/ThermalDynamics.cpp"
#endif

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

//NEEDS TO BE FIRST BECAUSE OF CALIBRATED COOLING ESTIMATES
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
#include "TOPOLOGY/IntegralEstimateMethod.cpp"

//////////////////
//MISCELLANEOUS //
//////////////////

#include "IO/SaveConfiguration.cpp"

#if INPUT_FLAG==YES_FLAG
#include "IO/LoadConfiguration.cpp"
#endif


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
    
    void Init(){
        
        ///////////////////////////////////
        //SET NCS MEASUREMENT PARAMETERS //
        ///////////////////////////////////
        
        
        // SET COOLING FREQUENCY //
        CoolingFrequency=20; ChernSimonsNumber::CoolingMethod::StandardCoolingDepth=5;
        
        Cooling::BlockFrequency=2*ChernSimonsNumber::CoolingMethod::StandardCoolingDepth;
        
        ChernSimonsNumber::CoolingMethod::StandardBlockingLevel=0;
        
        CalibFrequency=5*CoolingFrequency;
        
        ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel=1;
        
        ChernSimonsNumber::CoolingMethod::CalibrationCoolingDepth=1000.0;
        
        Cooling::NCsOutputFrequency=INT(1.0/SQR(Qs)+0.5);
        
        
        // CONSISTENCY CHECK //
        if(ChernSimonsNumber::CoolingMethod::StandardBlockingLevel*Cooling::BlockFrequency*SQR(Qs*GradientFlow::SqrtDcTauOverSpacing)>ChernSimonsNumber::CoolingMethod::StandardCoolingDepth){
            
            std::cerr << "#ERROR -- INCONSISTENT COOLING PARAMETERS" << std::endl;
            exit(0);
            
        }
        
        // SET GAUGE FIXING AND SPECTRA MEASUREMENT FREQUENCY //
        GaugeFixingFrequency=40; GaugeFixingTotalSteps=1200;
        
        // SAVE CONFIGURATIONS
        ConfigSaveFrequency=INT(500.0/(Qs*Dynamics::dTau));
        
        /////////////////////
        // INITIALIZATIONS //
        /////////////////////
        
        //INITIALIZE LATTICE
        Lattice::Init();
        
        //INITIALIZE GAUGE TRANSFORMATIONS
        GaugeTransformation::Init();
        
        //INITIALIZE GAUGE TRANSFORMED FIELDS
        GaugeFixedVariables::Init();

        //INITIALIZE FOURIER SPACE
        FourierSpace::Init();
        
        //INITIALIZE COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Init();
        
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
        OutStream << "#Qs= " << Qs << std::endl;
        OutStream << "#n0= " << n0 << std::endl;
        OutStream << "#LATTICE DATA" << std::endl;
        OutStream << "#Nx,Ny,Nz= " << GLinks::U->N[0] << " " << GLinks::U->N[1] << " " << GLinks::U->N[2] << std::endl;
        OutStream << "#ax,ay,az= " << GLinks::U->a[0] << " " << GLinks::U->a[1] << " " << GLinks::U->a[2] << std::endl;
        OutStream << "#PARAMTER LIST" << std::endl;
        OutStream << "#dTau= " << Dynamics::dTau << std::endl;
        OutStream << "#COOLING PARAMETERS" << std::endl;
        OutStream << "#SqrtDcTauOverSpacing = " << GradientFlow::SqrtDcTauOverSpacing << std::endl;
        OutStream << "#CoolingFrequency = " << CoolingFrequency << std::endl;
        OutStream << "#StandardBlockingLevel= " << ChernSimonsNumber::CoolingMethod::StandardBlockingLevel << std::endl;
        OutStream << "#StandardCoolingDepth= " << ChernSimonsNumber::CoolingMethod::StandardCoolingDepth << std::endl;
        OutStream << "#BlockFrequency= " << Cooling::BlockFrequency << std::endl;
        OutStream << "#CALIBRATION PARAMETERS" << std::endl;
        OutStream << "#CalibFrequency= " << CalibFrequency << std::endl;
        OutStream << "#CalibrationCoolingDepth= " << ChernSimonsNumber::CoolingMethod::CalibrationCoolingDepth
        << std::endl;
        OutStream << "#CalibBlockingLevel= " << ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel << std::endl;
        OutStream << "#GAUGE FIXING PARAMETERS" << std::endl;
        OutStream << "#GaugeFixingFrequency = " << GaugeFixingFrequency << std::endl;
        OutStream << "#GaugeFixingTotalSteps= " << GaugeFixingTotalSteps << std::endl;
        #if IC_FLAG==THERMAL_FLAG
        OutStream << "#ThermalDynamics::beta= " << ThermalDynamics::beta << std::endl;
        #endif

        OutStream.close();
        
    }
    
    ////////////////////////////////////
    // CLASSICAL YANG-MILLS EVOLUTION //
    ////////////////////////////////////
    
    void Evolve(DOUBLE MaxTime){
        
        
        ///////////////////////
        // OUTPUT MANAGEMENT //
        ///////////////////////
        
        CreateInfoFile();
        
        std::ofstream EnergyOutStream,CoolNCsOutStream;
        
        #if NAIVENCS_FLAG==YES_FLAG
        std::ofstream NaiveNCsOutStream;
        NaiveNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"NaiveNCsDotID",MY_MPI_RNG_SEED,".txt").c_str());
        #endif
        
        EnergyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"EnergyID",MY_MPI_RNG_SEED,".txt").c_str());
        
        CoolNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CoolNCsID",MY_MPI_RNG_SEED,".txt").c_str());
        
        ///////////////////////////////
        // INITIAL STATE OBSERVABLES //
        ///////////////////////////////
        
        // SAVE INITIAL CONFIGURATION //
        IO::SaveConfiguration(StringManipulation::StringCast("UOutT",Dynamics::Time()).c_str(),StringManipulation::StringCast("EOutT",Dynamics::Time()).c_str());
        
        // MEASURE BULK OBSERVABLES //
        Observables::Bulk::Update();
        // MEASURE HARD SCALES //
        Observables::HardScales::Update();
        
        EnergyOutStream << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
        
        // RESET COULOMB GAUGE FIXING ALGORITHMS //
        CoulombGaugeFixing::Reset();
        
        // PERFORM COULOMB GAUGE FIXING //
        CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",Dynamics::Time()).c_str(),GaugeFixingTotalSteps);
        
        // MEASURE SPECTRA
        Observables::Spectra::Compute(StringManipulation::StringCast("SpectraT",Dynamics::Time()).c_str());
        
        // SET COULOMB GAUGE FIELDS AS DYNAMICAL VARIABLES //
        GaugeTransformation::Operations::SaveFields();
        
        // RESET CHERN SIMONS- MEASUREMENST//
        ChernSimonsNumber::DeltaNCsRealTime=DOUBLE(0.0); ChernSimonsNumber::DeltaNCsCoolRealTime=DOUBLE(0.0);
        ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0); ChernSimonsNumber::DeltaNCsPreviousCooling=DOUBLE(0.0);

        // START COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Start();
        
        // CREATE INITIAL OUTPUT //
        CoolNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
        
        // CALIBRATE //
        //ChernSimonsNumber::CoolingMethod::Calibrate();
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        while(Dynamics::Time()<MaxTime){
            
            // UPDATE GAUGE LINKS AND ELECTRIC FIELD VARIABLES //
            Dynamics::Update();
            
            // NAIVE NCS DOT
            #if NAIVENCS_FLAG==YES_FLAG
            NaiveNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsRealTime << std::endl;
            #endif
            
            // CHECK BULK OBSERVABLES //
            if(Dynamics::tSteps%20==0){
                
                #define ENERGY_MEASUREMENT 1
                
                // MEASURE BULK OBSERVABLES //
                Observables::Bulk::Update();
                // MEASURE HARD SCALES //
                Observables::HardScales::Update();
                
                EnergyOutStream << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
                
            }
            
            // COOLED CHERN SIMONS DERIVATIVE //
            if(Dynamics::tSteps%CoolingFrequency==0){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Update();
                
                // MEASURE DIFFERNCE IN CHERN SIMONS NUMBER  //
                CoolNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
            }
            
             // COOLING CALIBRATION //
             if(Dynamics::tSteps%CalibFrequency==0){
             
                 // UPDATE COOLING //
                 //ChernSimonsNumber::CoolingMethod::Calibrate();
             
             }
            /*
             // SAVE U AND E
             if(Dynamics::tSteps%ConfigSaveFrequency==0){
             // OUTPUT CONFIGURATIONS
             IO::SaveConfiguration(StringManipulation::StringCast("UOutT",Dynamics::Time()).c_str(),StringManipulation::StringCast("EOutT",Dynamics::Time()).c_str());
             }
             */
            
            // PERFORM GAUGE FIXING AND SET COULOMB GAUGE FIELDS //
            if(Dynamics::tSteps%GaugeFixingFrequency==0){
                
                // RESET COULOMB GAUGE FIXING ALGORITHMS //
                CoulombGaugeFixing::Reset();
             
                // PERFORM COULOMB GAUGE FIXING //
                CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",Dynamics::Time()).c_str(),GaugeFixingTotalSteps);
             
                // MEASURE SPECTRA
                Observables::Spectra::Compute(StringManipulation::StringCast("SpectraT",Dynamics::Time()).c_str());
             
                // SET COULOMB GAUGE FIELDS AS DYNAMICAL VARIABLES //
                GaugeTransformation::Operations::SaveFields();
                
                // COOL GAUGE FIXED FIELDS FOR NEXT COOLING UPDATE //
                
            }
            
            
            // SANITY CHECKS //
            if(Dynamics::tSteps%1000==0){
                
                //CHECK GAUSS LAW VIOLATION
                Observables::GaussLaw::CheckViolation();
                
                //CHECK UNITARITY VIOLATION
                Observables::Unitarity::CheckViolation();
                
            }
            
            
            
        }
        
        // CLOSE OUTPUT STREAMS //
        EnergyOutStream.close();    CoolNCsOutStream.close();
        #if NAIVENCS_FLAG==YES_FLAG
        NaiveNCsOutStream.close();
        #endif
        
    }
    
    
    //////////////////////////
    //SIMULATION PROCDEDURE //
    //////////////////////////
    
    void Run(INT MPI_RNG_SEED){
        
        ///////////
        // SETUP //
        ///////////
        
        //SET SEED //
        MY_MPI_RNG_SEED=1443154311;//MPI_RNG_SEED;//1443037688;//1442551674;//1443037688;//
        
        //INITIALIZE RANDOM NUMBER GENERATOR //
        RandomNumberGenerator::Init(MY_MPI_RNG_SEED);
        
        //INITIALIZE DYNAMICS //
        Dynamics::Reset();
        
        ////////////////////////////
        // SET INITIAL CONDITIONS //
        ////////////////////////////
        #if IC_FLAG==QP_FLAG
        InitialConditions::SetQuasiParticles();
        //RestartInitialConditions::SetQuasiParticles();
        
        //IO::LoadConfiguration(GLinks::U,EFields::E);
        //IO::LoadLinks(GLinks::U);
        #endif
        
        #if IC_FLAG==THERMAL_FLAG
        InitialConditions::SetZero();
        // THERMALIZE
        ThermalDynamics::Thermalize(20,50.0);
        #endif
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        std::cerr << "### STARTING TIME EVOLUTION ###" << std::endl;
        Evolve(1000.0);
         
    }
    
    
}

#ifndef ENERGY_MEASUREMENT
NO ENERGY MEASUREMENT
#endif


