//INCLUDE BASIC DEFINITIONS AND COMPILER FLAGS
#include "Definitions.cpp"

//SET GAUGE GROUP
#define SU_Nc_FLAG SU2_FLAG

//INCLUDE SIMULATION PROCEDURE
#include "CoolSphaleronSimulation.cpp" // FOR RUNS

// INCLUDE MPI SAMPLING //
#include "MPI/Basic.cpp"

// INCLUDE COMMANDLINE PROCESSING //
#include "IO/COMMANDLINE/cfile.c"

//MAIN
int main(int argc,char **argv){
    
    ///////////////////////////////
    //INITIALIZE MPI ENVIRONMENT //
    ///////////////////////////////
    
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //SET MPI ID's AND NUMBER OF NODES
    MPI_Comm_rank(MPI_COMM_WORLD,&MPIBasic::ID);
    MPI_Comm_size(MPI_COMM_WORLD,&MPIBasic::NumberOfNodes);
    
    
    ////////////////////////////////////
    // INITIALIZE OPEN MP ENVIRONMENT //
    ////////////////////////////////////
    
    //CHECK THREAD COUNT
    std::cerr << "#NUMBER OF THREADS " << omp_get_max_threads() << std::endl;
    
    //INITIALIZE THREADED FFTW
    int FFTW3_THREAD_STATUS=fftw_init_threads();
    
    std::cerr << "#FFTW THREAD STATUS " << FFTW3_THREAD_STATUS << std::endl;
    
    if(FFTW3_THREAD_STATUS==1){
        
        fftw_plan_with_nthreads(omp_get_max_threads());
        
    }
    
    //////////////////////////////////
    //PROCESS COMMANDLINE ARGUMENTS //
    //////////////////////////////////
    
    INT NumberOfConfigurations=1;
    
    //////////////////////////////////
    //PROCESS COMMANDLINE ARGUMENTS
    //////////////////////////////////
    
    Konfig arguments(argc,argv);
    
    //GET OUTPUT FOLDER
    arguments.Getval("nconfs",NumberOfConfigurations);
    
    //////////////////////////
    // SET OUTPUT DIRECTORY //
    //////////////////////////
    
    char OutDir[256]="OUTPUT";
    
    arguments.Getval("o",OutDir);
    
    IO::SetOutputDirectory(OutDir);
    
#if INPUT_FLAG==YES_FLAG
    
    /////////////////////////
    // SET INPUT DIRECTORY //
    /////////////////////////
    
    char InDir[256]="INPUT";
    
    arguments.Getval("i",InDir);
    
    IO::SetInputDirectory(InDir);
    
    /////////////////////
    // SET INPUT FILES //
    /////////////////////
    
    // FOR LOADING FILES
    INT InputFileTime=0;
    INT InputFileID=1438889926;
    
    arguments.Getval("iT",InputFileTime);
    arguments.Getval("iID",InputFileID);
    
    IO::SetInputFile(InputFileTime,InputFileID);
    
#endif
    
    
    ////////////////////////////
    // DETERMINE LATTICE SIZE //
    ////////////////////////////
    
    INT NSites=-1;
    
    arguments.Getval("N",NSites);
    
    if(NSites>0){
        
        Lattice::N[0]=NSites;
        Lattice::N[1]=NSites;
        Lattice::N[2]=NSites;
        
        Lattice::Volume=NSites*NSites*NSites;
        
        std::cerr << "## LATTICE SIZE IS " << Lattice::N[0] << "x" << Lattice::N[1] << "x" << Lattice::N[2] << std::endl;
        
    }
    else{
        std::cerr << "## NUMBER OF SITES NOT SPECIFIED -- USING " << Lattice::N[0] << "x" << Lattice::N[1] << "x" << Lattice::N[2] << std::endl;
    }
    
    ///////////////////////////////
    // GET SIMULATION PARAMETERS //
    ///////////////////////////////
    
    arguments.Getval("Qs",Qs);
    arguments.Getval("n0",n0);
    
    std::cerr << "#Qs=" << Qs << " n0=" << n0 << std::endl;
    
    //////////////
    // SIMULATE //
    //////////////
    
    
    //COMMAND LINE OUTPUT
    std::cerr << "#GAUGE GROUP IS SU(" << Nc << ")" << std::endl;
    
    std::cerr << "#PRECISION IS " << MAX_DIGITS_PRECISION << " DIGITS" << std::endl;
    
    //INITIALIZE SIMULATION
    Simulation::Init();
    
    // SAMPLE DIFFERENT CONFIGURATIONS //
    for(INT n=0;n<NumberOfConfigurations;n++){
        
        //SET GLOBAL RANDOM NUMBER SEED//
        INT GLOBAL_RNG_SEED;
        
        if(MPIBasic::ID==0){
            GLOBAL_RNG_SEED=time(0);
        }
        
        // BROADCAST GLOBAL RANDOM SEED //
        MPI_Bcast(&GLOBAL_RNG_SEED, 1, MPI_INT,0,MPI_COMM_WORLD);
        
        // PERFORM CLASSICAL STATISTICAL SIMULATION //
        Simulation::Run(GLOBAL_RNG_SEED+MPIBasic::ID);
        
        // COMMADNLINE NOTIFICATION //
        std::cerr << "#COMPLETED " << GLOBAL_RNG_SEED+MPIBasic::ID << std::endl;
        
    }
    
    
    //SYNCHRONIZE ALL MPI NODES
    MPI_Barrier(MPI_COMM_WORLD);
    
    //FINALIZE MPI
    MPI_Finalize();
    
    //EXIT
    exit(0);
    
}
