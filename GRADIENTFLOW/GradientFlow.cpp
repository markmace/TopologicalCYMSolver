#ifndef __GRADIENTFLOW__CPP__
#define __GRADIENTFLOW__CPP__

namespace GradientFlow{
    
    
    //////////////////////
    //   COOLING TIME   //
    //////////////////////
    
    DOUBLE cTau;
    
    INT CoolingSteps;
    
    DOUBLE SqrtDcTauOverSpacing=0.5/D_SQRT2;
    
    
    // RESET //
    void Reset(){
        
        // RESET COOLING TIME AND STEPS //
        cTau=0.0; CoolingSteps=0;
        
    }
    
    // COOLING TIME //
    DOUBLE CoolingTime(){
        return cTau;
    }
    
    /////////////////////////////////////////////////////////
    //   COMPUTE GRADIENT FLOW UPDATE OF THE GAUGE LINKS   //
    /////////////////////////////////////////////////////////
    
    void UpdateElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        // SET DYNAMICAL CONSTANTS -2d\tau \sqrt{-g} a^{3} g^{\mu\alpha} g^{\nu\alpha}/(a_{\mu}^2 a_{\nu}^2)
        DOUBLE cE[Lattice::Dimension];
        
        DOUBLE gamma=(-2.0)*Dynamics::MetricDeterminant*U->aCube*Lattice::aScale;
        
        for(INT mu=0;mu<Lattice::Dimension;mu++){
            cE[mu]=Dynamics::gUpMetric[mu]/SQR(E->a[mu]);
        }
        
        #pragma omp parallel
        {
            
            //ALLOCATE BUFFERS TO COMPUTE PLAQUETTES
            SET_ELEMENTARY_PLAQUETTE_BUFFERS();
            SET_NEIGHBORING_PLAQUETTE_BUFFERS();
            
            //ALLOCATE BUFFERS TO COMPUTE TRACES OF PLAQUETTES
            SET_ELEMENTARY_COLOR_TRACE_BUFFERS();
            SET_NEIGHBORING_COLOR_TRACE_BUFFERS();
            
            //COMPUTE THE UPDATE FOR ALL GAUGE LINKS
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
                        
                        //////////////////////////////////////////////////////////////////////////////
                        //GET GRADIENT FLOW LINK UPDATE \partial_{\tau} U= - \partial H/ \partial U //
                        //////////////////////////////////////////////////////////////////////////////
                        
                        // COMPUTE E=\partial_{tau} A_{i}=F_{tau i} WHICH IS THE SAME AS (REAL) TIME DERIVATIVE OF ELECTRIC FIELD//
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            E->Get(x,y,z,0,a)[0]=gamma*cE[0]*(cE[1]*(ReTrITaUxy[a]-ReTrITaUxMy[a])-cE[2]*(ReTrITaUzx[a]-ReTrITaUMzx[a]));
                            E->Get(x,y,z,1,a)[0]=gamma*cE[1]*(cE[2]*(ReTrITaUyz[a]-ReTrITaUyMz[a])-cE[0]*(ReTrITaUxy[a]-ReTrITaUMxy[a]));
                            E->Get(x,y,z,2,a)[0]=gamma*cE[2]*(cE[0]*(ReTrITaUzx[a]-ReTrITaUzMx[a])-cE[1]*(ReTrITaUyz[a]-ReTrITaUMyz[a]));
                            
                        }
                        
                    }
                }
            }
            
        } // END PARALLEL
        
    }
    
    ///////////////////////////////////
    //   PERFROM GAUGE LINK UPDATE   //
    ///////////////////////////////////
    
    DOUBLE MaxStepSize;
    
    void UpdateGaugeLinks(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        
        #pragma omp parallel
        {
            
            //SU(Nc) MATRIX BUFFER //
            SU_Nc_FUNDAMENTAL_FORMAT OldLink[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT LinkUpdate[SUNcGroup::MatrixSize];
            
            //SET CONSTANTS CONSTANTS//
            DOUBLE cU[Lattice::Dimension]; DOUBLE DcTau=SQR(SqrtDcTauOverSpacing)*std::pow(U->aCube,2.0/3.0);
            
            for(INT mu=0;mu<Lattice::Dimension;mu++){
                cU[mu]=(-1.0)*(Dynamics::gDownMetric[mu]*SQR(E->a[mu])*DcTau)/(Dynamics::MetricDeterminant*U->aCube*Lattice::aScale);
            }
        
            
            // UPDATE ALL GAUGE LINKS //
            #pragma omp for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(int mu=0;mu<Lattice::Dimension;mu++){
                            
                            // COMPUTE MATRIX EXPONENTIALS //
                            SUNcAlgebra::Operations::MatrixIExp(cU[mu],E->Get(x,y,z,mu,0),LinkUpdate);
                            
                            // COPY THE OLD LINK TO LOCAL BUFFER //
                            COPY_SUNcMatrix(OldLink,U->Get(x,y,z,mu));
                            
                            // COMPUTE THE UPDATED LINK //
                            SUNcGroup::Operations::UU(LinkUpdate,OldLink,U->Get(x,y,z,mu));
                        }
                        
                    }
                }
            }

            
        } // END PARALLEL
        
        
    }
    
    
    ////////////////////////
    //COMPLETE UPDATE STEP//
    ////////////////////////
    
    // STANDARD COOLING //
    void StandardUpdate(GaugeLinks *U,ElectricFields *E){
        
        DOUBLE DcTau=SQR(SqrtDcTauOverSpacing)*std::pow(U->aCube,2.0/3.0)/SQR(Lattice::aScale);
        
        Dynamics::SetMetric();
        
        UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E);
        
        UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
        CoolingSteps++;
        
        cTau+=DcTau;
    
        
    }
    
    // COOLING WITH DELTA NCS MEASUREMENT //
    void CalibrationUpdate(GaugeLinks *U,ElectricFields *E){
        
        DOUBLE DcTau=SQR(SqrtDcTauOverSpacing)*std::pow(U->aCube,2.0/3.0)/SQR(Lattice::aScale);
        
        Dynamics::SetMetric();
        
        UpdateElectricFields(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E);
        
        ChernSimonsNumber::DeltaNCsCooling+=DOUBLE(0.5)*DcTau*ChernSimonsNumber::NCsDot(U,E);
        
        UpdateGaugeLinks(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
        ChernSimonsNumber::DeltaNCsCooling+=DOUBLE(0.5)*DcTau*ChernSimonsNumber::NCsDot(U,E);
        
        CoolingSteps++;
        
        cTau+=DcTau;
        
       
        
    }
  
}

#endif

