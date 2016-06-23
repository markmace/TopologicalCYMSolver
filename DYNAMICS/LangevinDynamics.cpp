#ifndef __LANGEVINDYNAMICS_CPP__
#define __LANGEVINDYNAMICS_CPP__

namespace LangevinVariables{
    
    ElectricFields *Xi;
    
    //INITIALIZE VARIABLE
    void Init(){
        
        Xi=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        
    }
    
}

namespace LangevinDynamics{
    
    /////////////////////
    //   TEMPERATURE   //
    /////////////////////
    
    DOUBLE beta=ThermalDynamics::beta; // SET TO VALUE FROM THERMAL ICs
    
    ////////////////////////////
    //   COLOR CONDUCTIVITY   //
    ////////////////////////////
    DOUBLE sigmac=1.0; // SCALES OUT
    
    //////////////
    //   TIME   //
    //////////////
    
    DOUBLE tau0L=0.005;
    
    //DISCRETIZED TIME AND NUMBER OF TIME STEPS
    DOUBLE tauL;
    INT tauLSteps=0;
    
    //TIME INCREMENT
    static const DOUBLE dTauL=0.001;
    
    //GET EVOLUTION TIME
    DOUBLE Time(){
        return (tauLSteps*dTauL);
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
        
        gUpMetric[0]=1.0; gUpMetric[1]=1.0; gUpMetric[2]=DOUBLE(1.0)/SQR(tauL);
        
        gDownMetric[0]=1.0; gDownMetric[1]=1.0; gDownMetric[2]=SQR(tauL);
        
        MetricDeterminant=tauL;
        
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
    
    
    //////////////
    //   INIT   //
    //////////////
    
    void Reset(){
        
        tauL=tau0L;   tauLSteps=0;
        
        SetMetric();
        
    }
    
    
    // CHANGE //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                         COMPUTE UPDATE OF THE LATTICE GAUGE LINKS                                                                    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                                                      //
    //   U_{mu}(x+dTauL)= Exp(-ig a_{mu}/a^3 -g_{\mu\nu} (E^{\nu} - \partial H/ \partial A_\mu d\tau/2) dTau / sqrt(-g)) U_{mu}(x)          //
    //                                                                                                                                      //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void UpdateGaugeLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        //SET CONSTANTS CONSTANTS//
        DOUBLE cU[Lattice::Dimension];
        
        for(INT mu=0;mu<Lattice::Dimension;mu++){
            cU[mu]=-dTauL*gDownMetric[mu]*SQR(U->a[mu])/(MetricDeterminant*U->aCube);
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
        
        //SET DYNAMICAL CONSTANTS -2d\tauL \sqrt{-g} a^{3} g^{\mu\alpha} g^{\nu\alpha}/(a_{\mu}^2 a_{\nu}^2)
        DOUBLE cE[Lattice::Dimension];
        DOUBLE cR[Lattice::Dimension];
        DOUBLE cS[Lattice::Dimension];
        
        DOUBLE gamma=-2.0*dTauL*MetricDeterminant*U->aCube;
        
        
        for(INT mu=0;mu<Lattice::Dimension;mu++){
            
            cE[mu]=gUpMetric[mu]/SQR(E->a[mu]);
            
            cR[mu]=gUpMetric[mu]*sqrt(dTauL)/SQR(LangevinVariables::Xi->a[mu]);

            cS[mu]=(gUpMetric[mu]*MetricDeterminant*U->aCube/SQR(E->a[mu]))*dTauL;

        }
        
        
        #pragma omp parallel
        {
            
            //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
            SET_ELEMENTARY_PLAQUETTE_BUFFERS();
            SET_NEIGHBORING_PLAQUETTE_BUFFERS();
            
            //ALLOCATE BUFFERS TO COMPUTE TRACES OF PLAQUETTES
            SET_ELEMENTARY_COLOR_TRACE_BUFFERS();
            SET_NEIGHBORING_COLOR_TRACE_BUFFERS();
            
            //BUFFERS FOR UPDATES
            DOUBLE EUpdate0,EUpdate1,EUpdate2;
            
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
                        // d/dt E^a_i(x,t)=-dH/dA^a_i-ColorConductivity*E^a_i(x,t)+Xi^a_i(x,t)
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            EUpdate0=gamma*cE[0]*(cE[1]*(ReTrITaUxy[a]-ReTrITaUxMy[a])-cE[2]*(ReTrITaUzx[a]-ReTrITaUMzx[a]))-cS[0]*E->Get(x,y,z,0,a)[0]+cR[0]*LangevinVariables::Xi->Get(x,y,z,0,a)[0];
                            
                            E->Get(x,y,z,0,a)[0]+=EUpdate0;
                            
                            EUpdate1=gamma*cE[1]*(cE[2]*(ReTrITaUyz[a]-ReTrITaUyMz[a])-cE[0]*(ReTrITaUxy[a]-ReTrITaUMxy[a]))-cS[1]*E->Get(x,y,z,1,a)[0]+cR[1]*LangevinVariables::Xi->Get(x,y,z,1,a)[0];
                            
                            E->Get(x,y,z,1,a)[0]+=EUpdate1;
                            
                            EUpdate2=gamma*cE[2]*(cE[0]*(ReTrITaUzx[a]-ReTrITaUzMx[a])-cE[1]*(ReTrITaUyz[a]-ReTrITaUMyz[a]))-cS[2]*E->Get(x,y,z,2,a)[0]+cR[2]*LangevinVariables::Xi->Get(x,y,z,2,a)[0];
                            
                            E->Get(x,y,z,2,a)[0]+=EUpdate2;
                        
                            /*
                            EUpdate0=-cS[0]*E->Get(x,y,z,0,a)[0]+cR[0]*LangevinVariables::Xi->Get(x,y,z,0,a)[0];
                            
                            E->Get(x,y,z,0,a)[0]+=EUpdate0;
                            
                            EUpdate1=-cS[1]*E->Get(x,y,z,1,a)[0]+cR[1]*LangevinVariables::Xi->Get(x,y,z,1,a)[0];
                            
                            E->Get(x,y,z,1,a)[0]+=EUpdate1;
                            
                            EUpdate2=-cS[2]*E->Get(x,y,z,2,a)[0]+cR[2]*LangevinVariables::Xi->Get(x,y,z,2,a)[0];
                            
                            E->Get(x,y,z,2,a)[0]+=EUpdate2;
                            */

                        }
                        
                    }
                }
            }
            
        } // END PARALLEL
    }
    
    // <Xi(x,t)^a_i Xi(y,t')^b_j>= 2 T \delta_ij \delta_xy \delta_ab \delta_tt' //
    void UpdateRandomGaussianNoise(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        DOUBLE SqrtTwoOverBeta=sqrt(2.0/beta);

        // SAMPLE GAUSSIAN ELECTRIC FIELDS
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    //SET ELECTRIC FIELDS TO GAUSSIAN RANDOM NUMBER DISTRIBUTION GIVEN AN INVERSE TEMPERATURE
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            LangevinVariables::Xi->Get(x,y,z,mu,a)[0]=RandomNumberGenerator::Gauss(SqrtTwoOverBeta);
                            
                        }
                    }
                    
                }
            }
        }
        
        // PROJECT GAUSS LAW //
        // OPTION TO RESTORE GAUSS LAW -- NOT USED IN IMPLEMENTATION FROM 1101.1167
        //GaussLawRestoration::Restore();
        // END OPTION
    }
    
    //OVERLOAD
    void UpdateRandomGaussianNoise(){
        
        UpdateRandomGaussianNoise(0,LangevinVariables::Xi->N[0]-1,0,LangevinVariables::Xi->N[1]-1,0,LangevinVariables::Xi->N[2]-1);
        
    }
    
    
    ////////////////////////
    //COMPLETE UPDATE STEP//
    ////////////////////////
    
    void UpdateFields(GaugeLinks *U,ElectricFields *E){
        
        SetMetric();
        
        UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E);
        
        UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
        tauLSteps++;
        
        
    }
    
    void UpdateFields(){
        
        UpdateFields(GLinks::U,EFields::E);
        
    }
    
    
    void Update(){
        
        // OUTPUT STREAM //
        std::ofstream EnergyOutStream;
        
        // DETERMINE RANDOM GAUSSIAN NOISE //
        UpdateRandomGaussianNoise();
        
        // EVOLVE FIELDS //
        UpdateFields();
                
        if(tauLSteps%500==0){
            
            //CHECK GAUSS LAW VIOLATION //
            Observables::GaussLaw::CheckViolation();
            
            //CHECK UNITARITY VIOLATION //
            Observables::Unitarity::CheckViolation();
            
            
        }
                
        // RESET TIME
        // OPTION TO RESET THE LANGEVIN TIME
        //Reset();
        // END OPTION
        
    }
    
    
}

#endif
