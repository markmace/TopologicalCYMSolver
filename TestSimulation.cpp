//////////////////////////////////////////
// CHERN SIMONS NUMBER CALCULATION TEST //
//////////////////////////////////////////

#define METRIC_FLAG MINKOWSKI_FLAG
#define NAIVENCS_FLAG NO_FLAG
#define INPUT_FLAG NO_FLAG

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

//INCLUDE MACROS TO COMPUTE GAUSS LAW VIOLATION
#include "CALC/GaussViolation.cpp"

//INCLUDE MACROS TO PERFORM GAUGE TRANSFORMATIONS
#include "CALC/GaugeTransformation.cpp"

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
#include "SCALEOBSERVABLES/HardScales.cpp"

//INCLUDE MEASUREMENT OF WILSON LOOPS
#include "SCALEOBSERVABLES/WilsonLoop.cpp"
#include "SCALEOBSERVABLES/TimeWilsonLoop.cpp"

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

////////////////////////////
//GAUGE FIXING PROCEDURES //
////////////////////////////

//INCLUDE GAUGE TRANSFORMATION RULES AND BUFFERS
#include "GAUGETRANSFORMATION/GaugeTransformation.cpp"

//INCLUDE COULOMB GAUGE FIXING PROCEDURE
#include "GAUGETRANSFORMATION/CoulombGaugeFixing.cpp"

//INCLUDE MEASUREMENT OF SPECTRA
#include "OBSERVABLES/Spectra.cpp"

//INCLUDE MEASUREMENT OF E^2,B^2,E.B^2 CORRELATORS
#include "OBSERVABLES/Correlators.cpp"

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
#include "TOPOLOGY/SlaveFieldMethod.cpp"
#include "TOPOLOGY/IntegralEstimateMethod.cpp"

///////////////////////
//INITIAL CONDITIONS //
///////////////////////

#include "INITIALCONDITIONS/SetPureGauge.cpp"

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


namespace Test {
    
    // MPI RANDOM NUMBER SEED //
    INT MY_MPI_RNG_SEED;
    
    /////////////////////////////
    //INITIALIZATION AND EXIT  //
    /////////////////////////////
    
    INT CoolingFrequency, GaugeFixingFrequency, GaugeFixingTotalSteps,ConfigSaveFrequency,SpectraSaveFrequency;
    
    INT RefTimeFrequency, MeasureFrequency;
    
    void Init(){
        
        ///////////////////////////////////
        //SET NCS MEASUREMENT PARAMETERS //
        ///////////////////////////////////
        
        
        // SET COOLING FREQUENCY //
        CoolingFrequency=1;//INT(0.1/(Qs*Dynamics::dTau));
        
        ChernSimonsNumber::CoolingMethod::StandardCoolingDepth=1.0/SQR(Qs);
        
        Cooling::BlockFrequency=20*CoolingFrequency;//INT(0.2/SQR(DOUBLE(Qs*GradientFlow::SqrtDcTauOverSpacing))+0.5);
        
        Cooling::NCsOutputFrequency=INT(1.0/SQR(Qs)+0.5);
        
        ChernSimonsNumber::CoolingMethod::StandardBlockingLevel=0;
        
        // CONSISTENCY CHECK //
        if(ChernSimonsNumber::CoolingMethod::StandardBlockingLevel*Cooling::BlockFrequency*SQR(Qs*GradientFlow::SqrtDcTauOverSpacing)>ChernSimonsNumber::CoolingMethod::StandardCoolingDepth){
            
            std::cerr << "#ERROR -- INCONSISTENT COOLING PARAMETERS" << std::endl;
            exit(0);
            
        }
        
        // SET GAUGE FIXING AND SPECTRA MEASUREMENT FREQUENCY //
        GaugeFixingFrequency=INT(1.0/(Qs*Dynamics::dTau));
        
        GaugeFixingTotalSteps=1000;
        
        // SPECTRA
        SpectraSaveFrequency=INT(100.0/(Qs*Dynamics::dTau));
        
        // SAVE CONFIGURATIONS
        ConfigSaveFrequency=INT(500.0/(Qs*Dynamics::dTau));
        ////////////////////////////////////////
        // SPACE-TIME WILSON LOOP MEASUREMENT //
        ////////////////////////////////////////
        
        // SPACE-TIME WILSON LOOP REFERENCE
        RefTimeFrequency=INT(100.0/(Qs*Dynamics::dTau));
        
        // SPACE-TIME WILSON LOOP MEASUREMENT
        MeasureFrequency=INT(0.5/(Qs*Dynamics::dTau));
        
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
        
        //INITIALIZE SLAVE FIELD METHOD
        ChernSimonsNumber::SlaveFieldMethod::Init();
        
        //INITIALIZE COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Init();
        
        //INITIALIZE FOR CORRELATOR MEASUREMENT //
        Observables::Correlators::Init();
        
        //INITIALIZE FOR WILSON LOOP MEASUREMENT //
        Observables::WilsonLoop::Init();
        
        //INITIALIZE FOR TIME WILSON LOOP MEASUREMENT //
        Observables::TimeWilsonLoop::Init();
        
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
        OutStream << "#SqrtDcTau = " << GradientFlow::SqrtDcTau << std::endl;
        OutStream << "#CoolingFrequency = " << CoolingFrequency << std::endl;
        OutStream << "#StandardBlockingLevel= " << ChernSimonsNumber::CoolingMethod::StandardBlockingLevel << std::endl;
        OutStream << "#StandardCoolingDepth= " << ChernSimonsNumber::CoolingMethod::StandardCoolingDepth << std::endl;
        OutStream << "#BlockFrequency= " << Cooling::BlockFrequency << std::endl;
        OutStream << "#GAUGE FIXING PARAMETERS" << std::endl;
        OutStream << "#GaugeFixingFrequency = " << GaugeFixingFrequency << std::endl;
        OutStream << "#GaugeFixingTotalSteps= " << GaugeFixingTotalSteps << std::endl;
        
        OutStream.close();
        
    }
    
    ////////////////////////////////////
    // CLASSICAL YANG-MILLS EVOLUTION //
    ////////////////////////////////////
    
    void Compute(){
        
        std::ofstream CoolNCsOutStream;
        
        CoolNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CoolNCsID",MY_MPI_RNG_SEED,".txt").c_str());
        
        ///////////////////////
        // OUTPUT MANAGEMENT //
        ///////////////////////
        
        CreateInfoFile();
        
        ///////////////////////////////////
        // LINK AND FIELD CONFIGURATIONS //
        ///////////////////////////////////
        
        GaugeLinks *U1=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        GaugeLinks *U2=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
            
        ElectricFields *E1=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        ElectricFields *E2=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        
        // GENERATE FIRST PURE GAUGE CONFIGURATION //
        INT BAM=0;
        
        while(BAM==0){
        
        
            std::cerr << "#SETTING FIRST SET OF CONFIGURATIONS " << std::endl;
            InitialConditions::SetPureGauge();
            
            BAM=ChernSimonsNumber::SlaveFieldMethod::NCS();
            
            std::cerr << BAM << std::endl;
            
        }
    
        // COMPUTE TOPOLOGY //
        DOUBLE SlaveNCS0=ChernSimonsNumber::SlaveFieldMethod::NCS();
        DOUBLE VacEstNCS0=ChernSimonsNumber::VacuumEstimator::NCS(GLinks::U);
        DOUBLE IntEstNCS0=ChernSimonsNumber::IntegralEstimateMethod::NCS();
        std::cerr << "#SlaveNCS0 " << SlaveNCS0 << " " << VacEstNCS0 << " " << IntEstNCS0 << std::endl;
        
        // RESET CHERN SIMONS- MEASUREMENST//
        ChernSimonsNumber::DeltaNCsRealTime=DOUBLE(0.0); ChernSimonsNumber::DeltaNCsCoolRealTime=DOUBLE(0.0);  ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
        
        // START COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Start();
        
        // CREATE INITIAL OUTPUT //
        std::cerr << "COOLNCS " << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
        
        // RESET COULOMB GAUGE FIXING ALGORITHMS //
        CoulombGaugeFixing::Reset();
        // PERFORM GAUGE FIXING //
        CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",Dynamics::Time()).c_str(),1000);
        
        // SAVE LINKS AND FIELDS //
        //Copy(U1,GLinks::U);
        //Copy(E1,EFields::E);
        
        // COMPUTE TOPOLOGY //
        DOUBLE SlaveNCS1=ChernSimonsNumber::SlaveFieldMethod::NCS();
        DOUBLE VacEstNCS1=ChernSimonsNumber::VacuumEstimator::NCS(GaugeFixedVariables::U);
        DOUBLE IntEstNCS1=ChernSimonsNumber::IntegralEstimateMethod::NCS();
        std::cerr << "#SlaveNCS1 " << SlaveNCS1 << " " << VacEstNCS1 << " " << IntEstNCS1 << std::endl;
        
        // GENERATE SECOND PURE GAUGE CONFIGURATION //
        std::cerr << "#SETTING SECOND SET OF CONFIGURATIONS " << std::endl;
        InitialConditions::SetPureGauge();
        
        // UPDATE COOLING //
        ChernSimonsNumber::CoolingMethod::Update();
        
        // MEASURE DIFFERNCE IN CHERN SIMONS NUMBER  //
        std::cerr << "COOLNCS " << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
        
        // SAVE LINKS AND FIELDS //
        //Copy(U2,GLinks::U);
        //Copy(E2,EFields::E);
        
        // COMPUTE TOPOLOGY //
        DOUBLE SlaveNCS2=ChernSimonsNumber::SlaveFieldMethod::NCS();
        DOUBLE VacEstNCS2=ChernSimonsNumber::VacuumEstimator::NCS(GLinks::U);
        DOUBLE IntEstNCS2=ChernSimonsNumber::IntegralEstimateMethod::NCS();
        std::cerr << "#SlaveNCS2 " << SlaveNCS2 << " " << VacEstNCS2 << " " << IntEstNCS2 << std::endl;
        
        // RESET COULOMB GAUGE FIXING ALGORITHMS //
        CoulombGaugeFixing::Reset();
        // PERFORM GAUGE FIXING //
        CoulombGaugeFixing::SetCoulombGauge(StringManipulation::StringCast("GaugeFixingT",Dynamics::Time()).c_str(),1000);
        
        // COMPUTE TOPOLOGY //
        DOUBLE SlaveNCS3=ChernSimonsNumber::SlaveFieldMethod::NCS();
        DOUBLE VacEstNCS3=ChernSimonsNumber::VacuumEstimator::NCS(GaugeFixedVariables::U);
        DOUBLE IntEstNCS3=ChernSimonsNumber::IntegralEstimateMethod::NCS();
        std::cerr << "#SlaveNCS3 " << SlaveNCS3 << " " << VacEstNCS3 << " " << IntEstNCS3 << std::endl;
        
        delete U1;
        delete U2;
        delete E1;
        delete E2;
        
        CoolNCsOutStream.close();
        
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
        
        //INITIALIZE DYNAMICS //
        Dynamics::Reset();
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        std::cerr << "### STARTING COMPUTATIONS ###" << std::endl;
        Compute();
        
        
    }
    
    
}

