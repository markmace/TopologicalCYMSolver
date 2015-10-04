///////////////////////////////////////////
//INCLUDE FREE FIELD POLARIZATION VECTORS//
///////////////////////////////////////////

#if METRIC_FLAG==BJORKEN_FLAG
#include "../MISC/BjorkenWaveVectors.cpp"
#endif

#if METRIC_FLAG==MINKOWSKI_FLAG
#include "../MISC/MinkowskiWaveVectors.cpp"
#endif

///////////////////////////////////////////
//    INCLUDE FOURIER SPACE VARIABLES    //
///////////////////////////////////////////

#include "../LATTICE/FourierSpace.cpp"

/////////////////////////////////////////////
//    QUASI PARTICLE INITIAL CONDITIONS    //
/////////////////////////////////////////////


namespace InitialConditions{
    
    // INITIAL SINGLE PARTICLE DISTRIBUTION //
    DOUBLE SingleParticleDistribution(DOUBLE pXValue,DOUBLE pYValue,DOUBLE pZValue){
        
        DOUBLE pSqr=SQR(pXValue)+SQR(pYValue)+SQR(pZValue);
       
        if(pSqr<SQR(Qs)){
            
            //return DOUBLE(n0);
            return n0/sqrt(pSqr/SQR(Qs)+0.1);
        }
        else{
            return 0.0;
        }
        
    }
    
    // SET FIELD VARIABLES IN MOMENTUM SPACE //
    void SetMomentumSpaceVariables(){
        
        // MOMENTUM VARIABLES //
        COMPLEX cDpx,cDpy,cDpz; DOUBLE pXValue,pYValue,pZValue;        
        // POLARIZATION VECTORS //
        COMPLEX Xi1[Lattice::Dimension]; COMPLEX Xi1Dot[Lattice::Dimension];
        COMPLEX Xi2[Lattice::Dimension]; COMPLEX Xi2Dot[Lattice::Dimension];
        
        // MODE OCCUPANCY //
        DOUBLE OccupancyFactor; DOUBLE NormalizationFactor=sqrt(GLinks::U->Volume*Lattice::aCube);
        
        // STOCHASTIC VARIABLES //
        COMPLEX c1,c2; //DOUBLE g1,g2,phi1,phi2;
        
        /////////////////////////////////////////////////////////////////
        // SET SUPERPOSITION OF QUASI PARTICLE MODES IN MOMENTUM SPACE //
        /////////////////////////////////////////////////////////////////
        
        for(INT pXIndex=0;pXIndex<=GLinks::U->N[0]-1;pXIndex++){
            for(INT pYIndex=0;pYIndex<=GLinks::U->N[1]-1;pYIndex++){
                for(INT pZIndex=0;pZIndex<=GLinks::U->N[2]-1;pZIndex++){
                    
                    // DETERMINE MOMENTA //
                    GetMomenta(pXIndex,pYIndex,pZIndex,cDpx,cDpy,cDpz);
                    
                    pXValue=SIGN(real(cDpx))*sqrt(Dynamics::gUpMetric[0]*SQR(Lattice::aScale)*SQR_ABS(cDpx));
                    pYValue=SIGN(real(cDpy))*sqrt(Dynamics::gUpMetric[1]*SQR(Lattice::aScale)*SQR_ABS(cDpy));
                    pZValue=SIGN(real(cDpz))*sqrt(Dynamics::gUpMetric[2]*SQR(Lattice::aScale)*SQR_ABS(cDpz));
                                        
                    // DETERMINE POLARIZATION VECTORS //
                    PolarizationVectors::Compute(pXIndex,cDpx,pYIndex,cDpy,pZIndex,cDpz,Xi1,Xi1Dot,Xi2,Xi2Dot);
                    
                    // DETERMINE MODE OCCUPANCY //
                    OccupancyFactor=sqrt(SingleParticleDistribution(pXValue,pYValue,pZValue));
                    
                    // SET MOMENTUM SPACE MODES //
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        // SAMPLE STOCHASTIC COEFFICIENTS //
                        /*g1Re=RandomNumberGenerator::Gauss();  phi1=DOUBLE(2.0)*PI*RandomNumberGenerator::rng();
                        g2Re=RandomNumberGenerator::Gauss();  phi2=DOUBLE(2.0)*PI*RandomNumberGenerator::rng();
                            
                        // SET COMPLEX GAUSSIAN RANDOM NUMBERS //
                        c1=g1*exp(ComplexI*phi1); c2=g2*exp(ComplexI*phi2);*/
                        
                        c1=RandomNumberGenerator::ComplexGauss(); c2=RandomNumberGenerator::ComplexGauss();
                        
                        // SET MODE OCCUPANCIES
                        FourierSpace::A0->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(c1*Xi1[0]+c2*Xi2[0]));
                        FourierSpace::A1->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(c1*Xi1[1]+c2*Xi2[1]));
                        FourierSpace::A2->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(c1*Xi1[2]+c2*Xi2[2]));
                        
                        FourierSpace::E0->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(Dynamics::MetricDeterminant*Dynamics::gUpMetric[0])*(c1*Xi1Dot[0]+c2*Xi2Dot[0]));
                        FourierSpace::E1->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(Dynamics::MetricDeterminant*Dynamics::gUpMetric[1])*(c1*Xi1Dot[1]+c2*Xi2Dot[1]));
                        FourierSpace::E2->SetP(pXIndex,pYIndex,pZIndex,a,OccupancyFactor*NormalizationFactor*(Dynamics::MetricDeterminant*Dynamics::gUpMetric[2])*(c1*Xi1Dot[2]+c2*Xi2Dot[2]));
                        
                        
                    }
                    
                }
            }
        }
        
    }
    
    ////////////////////////////////////
    // PERFORM FAST FOURIER TRANSFORM //
    ////////////////////////////////////
    
    void PerformFourierTransform(){
        
        FourierSpace::A0->ExecutePtoX();  FourierSpace::A1->ExecutePtoX();  FourierSpace::A2->ExecutePtoX();
        FourierSpace::E0->ExecutePtoX();  FourierSpace::E1->ExecutePtoX();  FourierSpace::E2->ExecutePtoX();
        
    }
    
    ////////////////////////////////////
    // SET COORDINATE SPACE VARIABLES //
    ////////////////////////////////////
    
    void SetCoordinateSpaceVariables(){
        
        // FOURIER SPACE NORMALIZATION FACTOR //
        DOUBLE NormalizationFactor=DOUBLE(1.0)/(GLinks::U->Volume*Lattice::aCube);
        
        // GAUGE FIELD BUFFERS //
        DOUBLE A0Buffer[SUNcAlgebra::VectorSize]; 
        DOUBLE A1Buffer[SUNcAlgebra::VectorSize];
        DOUBLE A2Buffer[SUNcAlgebra::VectorSize];
        
        for(INT x=0;x<=GLinks::U->N[0]-1;x++){
            for(INT y=0;y<=GLinks::U->N[1]-1;y++){
                for(INT z=0;z<=GLinks::U->N[2]-1;z++){
                    
                    // SET GAUGE LINK VARIABLES //
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        A0Buffer[a]=NormalizationFactor*GLinks::U->a[0]*DOUBLE(2.0)*real(FourierSpace::A0->GetX(x,y,z,a));
                        A1Buffer[a]=NormalizationFactor*GLinks::U->a[1]*DOUBLE(2.0)*real(FourierSpace::A1->GetX(x,y,z,a));
                        A2Buffer[a]=NormalizationFactor*GLinks::U->a[2]*DOUBLE(2.0)*real(FourierSpace::A2->GetX(x,y,z,a));
                        
                    }
                    
                    SUNcAlgebra::Operations::MatrixIExp(-1.0,A0Buffer,GLinks::U->Get(x,y,z,0));
                    SUNcAlgebra::Operations::MatrixIExp(-1.0,A1Buffer,GLinks::U->Get(x,y,z,1));
                    SUNcAlgebra::Operations::MatrixIExp(-1.0,A2Buffer,GLinks::U->Get(x,y,z,2));
                    
                    // SET ELECTRIC FIELD VARIABLES //
                    for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                        
                        EFields::E->Get(x,y,z,0,a)[0]=NormalizationFactor*(Lattice::aCube/EFields::E->a[0])*DOUBLE(2.0)*real(FourierSpace::E0->GetX(x,y,z,a));
                        EFields::E->Get(x,y,z,1,a)[0]=NormalizationFactor*(Lattice::aCube/EFields::E->a[1])*DOUBLE(2.0)*real(FourierSpace::E1->GetX(x,y,z,a));
                        EFields::E->Get(x,y,z,2,a)[0]=NormalizationFactor*(Lattice::aCube/EFields::E->a[2])*DOUBLE(2.0)*real(FourierSpace::E2->GetX(x,y,z,a));

                    }
                                        
                }
            }
        }

        
    }
    
    
    
    void SetQuasiParticles(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "#SETTING QUASI-PARTICLE INITIAL CONDITIONS" << std::endl;
        
        SetMomentumSpaceVariables();
        
        PerformFourierTransform();
        
        SetCoordinateSpaceVariables();

        ///////////////////////////////
        // CHECK INITIAL OBSERVABLES //
        ///////////////////////////////
        
        std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
        
        //CHECK GAUSS LAW VIOLATION //
        Observables::GaussLaw::CheckViolation();
        
        //CHECK UNITARITY VIOLATION //
        Observables::Unitarity::CheckViolation();
        
        //CHECK INITIAL ENERGY DENSITY //
        Observables::Bulk::Update();

        std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
        std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
        
        // RESTORE GAUSS LAW //
        GaussLawRestoration::Restore();
        
        std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
        
        //CHECK GAUSS LAW VIOLATION //
        Observables::GaussLaw::CheckViolation();
        
        //CHECK UNITARITY VIOLATION //
        Observables::Unitarity::CheckViolation();
        
        //CHECK INITIAL ENERGY DENSITY //
        Observables::Bulk::Update();
        
        std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
        std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
            
    }
    
    
    
}