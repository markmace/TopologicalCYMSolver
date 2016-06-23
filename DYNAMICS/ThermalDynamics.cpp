#ifndef __THERMALDYNAMICS_CPP__
#define __THERMALDYNAMICS_CPP__

namespace ThermalDynamics{
    
    /////////////////////
    //   TEMPERATURE   //
    /////////////////////
    
    DOUBLE beta=2.0; // --> MULTIPLE BY 2*Nc TO GET GDM VALUE
    
    //////////////
    //   TIME   //
    //////////////
    
    DOUBLE tau0=0.005;
    
    //DISCRETIZED TIME AND NUMBER OF TIME STEPS
    DOUBLE tau;
    INT tSteps=0;
    
    //TIME INCREMENT
    static const DOUBLE dTau=0.0375;
    
    //GET EVOLUTION TIME
    DOUBLE Time(){
        return (tSteps*dTau);
    }
    
    ////////////////
    //   METRIC   //
    ////////////////
    
    
    //DIAGONAL SPATIAL COMPONENTS OF THE METRIC -g^{mu\nu}
    DOUBLE gUpMetric[Lattice::Dimension];
    
    //DIAGONAL SPATIAL COMPONENTS OF THE METRIC -g_{mu\nu}
    DOUBLE gDownMetric[Lattice::Dimension];
    
    //METRIC DETERMINANT sqrt{-g_{\mu\nu}(x)}
    DOUBLE MetricDeterminant;
    
    
    //SET MINKOWSKI METRIC
    void SetMinkowskiMetric(){
        
        gUpMetric[0]=1.0; gUpMetric[1]=1.0; gUpMetric[2]=1.0;
        
        gDownMetric[0]=1.0; gDownMetric[1]=1.0; gDownMetric[2]=1.0;
        
        MetricDeterminant=1.0;
        
    }
    
    //SET BJORKEN METRIC
    void SetBjorkenMetric(){
        
        gUpMetric[0]=1.0; gUpMetric[1]=1.0; gUpMetric[2]=DOUBLE(1.0)/SQR(tau);
        
        gDownMetric[0]=1.0; gDownMetric[1]=1.0; gDownMetric[2]=SQR(tau);
        
        MetricDeterminant=tau;
                
    }
    
    //SET METRIC
    void SetMetric(){
        
        #if METRIC_FLAG==MINKOWSKI_FLAG
        SetMinkowskiMetric();
        #endif
        
        #if METRIC_FLAG==BJORKEN_FLAG
        SetBjorkenMetric();  
        #endif
        
    }
    
    
    ////////////////
    //   INIT   //
    ////////////////
    
    void Reset(){
        
        tau=tau0;   tSteps=0;
        
        SetMetric();
        
    }

    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                         COMPUTE UPDATE OF THE LATTICE GAUGE LINKS                            //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    //   U_{mu}(x+dTau)= Exp(-ig a_{mu}/a^3 -g_{\mu\nu} E^{\nu} dTau / sqrt(-g)) U_{mu}(x)          //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
 
    void UpdateGaugeLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        //SET CONSTANTS CONSTANTS//
        DOUBLE cU[Lattice::Dimension];
        
        for(INT mu=0;mu<Lattice::Dimension;mu++){
            cU[mu]=-dTau*gDownMetric[mu]*SQR(U->a[mu])/(MetricDeterminant*U->aCube);
        }
        
        #pragma omp parallel
        {
            
            //MATRIX BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT LinkUpdate[SUNcGroup::MatrixSize];       
            SU_Nc_FUNDAMENTAL_FORMAT OldLink[SUNcGroup::MatrixSize];
            
            //UPDATE ALL GAUGE LINKS
            #pragma omp for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //UPDATE ALL LORENTZ COMPONENTS
                        for(INT mu=0;mu<Lattice::Dimension;mu++){
                            
                            //COMPUTE MATRIX EXPONENTIAL
                            SUNcAlgebra::Operations::MatrixIExp(cU[mu],E->Get(x,y,z,mu,0),LinkUpdate);
                            
                            //COPY THE OLD LINK
                            COPY_SUNcMatrix(OldLink,U->Get(x,y,z,mu));
                            
                            //COMPUTE UPDATED LINK
                            SUNcGroup::Operations::UU(LinkUpdate,OldLink,U->Get(x,y,z,mu));
                        }
                        
                    }
                }
            }
            
        } // END PARALLEL
        
    }
    
    
    //////////////////////////////////////////////////////////
    //      COMPUTE UPDATE OF THE LATTICE GAUGE LINKS       //
    //////////////////////////////////////////////////////////
    
    void UpdateElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        //SET DYNAMICAL CONSTANTS -2d\tau \sqrt{-g} a^{3} g^{\mu\alpha} g^{\nu\alpha}/(a_{\mu}^2 a_{\nu}^2)
        DOUBLE cE[Lattice::Dimension];
        
        DOUBLE gamma=-2.0*dTau*MetricDeterminant*U->aCube;
        
        for(INT mu=0;mu<Lattice::Dimension;mu++){
            cE[mu]=gUpMetric[mu]/SQR(E->a[mu]);
        }
        
   
        #pragma omp parallel
        {
            
            //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
            SET_ELEMENTARY_PLAQUETTE_BUFFERS();
            SET_NEIGHBORING_PLAQUETTE_BUFFERS();
            
            //ALLOCATE BUFFERS TO COMPUTE TRACES OF PLAQUETTES
            SET_ELEMENTARY_COLOR_TRACE_BUFFERS();
            SET_NEIGHBORING_COLOR_TRACE_BUFFERS();
            
            //UPDATE ALL GAUGE LINKS
            #pragma omp for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //COMPUTE ELEMENTARY PLAQUETTES AND COLOR TRACES
                        COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z);
                        COMPUTE_ELEMENTARY_COLOR_TRACES();
                        
                        //COMPUTE NEIGHBORING PLAQUETTES AND COLOR TRACES
                        COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z);
                        COMPUTE_NEIGHBORING_COLOR_TRACES();
                        
                        
                        //UPDATE ELECTRIC FIELD VARIABLES
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            E->Get(x,y,z,0,a)[0]+=gamma*cE[0]*(cE[1]*(ReTrITaUxy[a]-ReTrITaUxMy[a])-cE[2]*(ReTrITaUzx[a]-ReTrITaUMzx[a]));
                            E->Get(x,y,z,1,a)[0]+=gamma*cE[1]*(cE[2]*(ReTrITaUyz[a]-ReTrITaUyMz[a])-cE[0]*(ReTrITaUxy[a]-ReTrITaUMxy[a]));
                            E->Get(x,y,z,2,a)[0]+=gamma*cE[2]*(cE[0]*(ReTrITaUzx[a]-ReTrITaUzMx[a])-cE[1]*(ReTrITaUyz[a]-ReTrITaUMyz[a]));
                            
                        }
                        
                    }
                }
            }
            
        } // END PARALLEL
    }
    
    void ReassignElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        // SAMPLE GAUSSIAN ELECTRIC FIELDS
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    //SET ELECTRIC FIELDS TO GAUSSIAN RANDOM NUMBER DISTRIBUTION GIVEN AN INVERSE TEMPERATURE
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            EFields::E->Get(x,y,z,mu,a)[0]=RandomNumberGenerator::Gauss()/sqrt(beta);
                        }
                    }
                    
                }
            }
        }
        
        // PROJECT GAUSS LAW //
        GaussLawRestoration::Restore();
        
    }
    
    //OVERLOAD
    void ReassignElectricFields(){
        
        ReassignElectricFields(0,EFields::E->N[0]-1,0,EFields::E->N[1]-1,0,EFields::E->N[2]-1);
    
    }
    

    ////////////////////////
    //COMPLETE UPDATE STEP//
    ////////////////////////
    
    void Update(GaugeLinks *U,ElectricFields *E){
        
        SetMetric();
        
        UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E);

        UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
        tSteps++;

        
    }
    
    void Update(){
        
        Update(GLinks::U,EFields::E);
        
    }
    
    
    void Thermalize(INT MaxSteps,DOUBLE MaxTime){
        
        std::ofstream EnergyOutStream;
    EnergyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"ThermalEnergy","ID",RandomNumberGenerator::MySEED,".txt").c_str());
    
        for(INT i=0;i<MaxSteps;i++){
            
            // COMMANDLINE READOUT //
            std::cerr << "#REASSIGNING ELECTRIC FIELDS " << std::endl;
            
            // REASSIGN ELECTRIC FIELDS TO THERMAL CONFIGURATION //
            ReassignElectricFields();
            
            // EVOLVE TO THERMALIZE GAUGE LINKS //
            while(Time()<MaxTime){
                
                // EVOLVE //
                Update();
                
                // MEASURE ENERGY //
                if(tSteps%20==0){
                    
                    Observables::Bulk::Update();
                
                    EnergyOutStream << Time()+i*MaxTime << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << std::endl;
                }
                
                if(tSteps%500==0){
                    Observables::Unitarity::CheckViolation();
                }
                
            }
            
            // RESET TIME
            Reset();
            
        }
        
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

#endif
