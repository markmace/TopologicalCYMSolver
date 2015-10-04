namespace InitialConditions{
    
    void SetPureGauge(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        // INITIALIZE GAUGE LINKS AND ELECTRIC FIELDS //
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                
                    // SET ELECTRIC FIELDS TO ZERO //
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            EFields::E->Get(x,y,z,mu,a)[0]=0.0;
                        }
                    }
                    
                    // SET GAUGE LINKS TO UNITY //
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        COPY_SUNcMatrix(GLinks::U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                    }

                    
                }
            }
        }
        
        // RANDOM PURE GAUGE TRANSFORM LINKS //
        GaugeTransformation::SetRandom();
        /*
        // CREATE OBJECTS TO BE TRANSFORMED //
        GaugeLinks *UGT=new GaugeLinks(GLinks::U->N[0],GLinks::U->N[1],GLinks::U->N[2],GLinks::U->a[0],GLinks::U->a[1],GLinks::U->a[2]);
        ElectricFields *EGT=new ElectricFields(EFields::E->N[0],EFields::E->N[1],EFields::E->N[2],EFields::E->a[0],EFields::E->a[1],EFields::E->a[2]);
        
        //std::cerr << "#TEST3 " << std::endl;
        
        // GAUGE TRANSFORMED OBJECTS //
        GaugeTransformation::Operations::GaugeTransformLinks(GLinks::U,UGT,GaugeTransformation::G);
        //GaugeTransformation::Operations::GaugeTransformElectricFields(EFields::E,EGT,GaugeTransformation::G);
        
        // COPY GAUGE TRANSFORMED OBJECTS
        Copy(GLinks::U,UGT);
        //Copy(EFields::E,EGT);
        
        // DELETE OLD OBJECTS //
        delete UGT;
        delete EGT;
        */

        GaugeTransformation::Operations::GaugeTransformLinks();
        GaugeTransformation::Operations::GaugeTransformElectricFields();
        
        Copy(GLinks::U,GaugeFixedVariables::U);
        Copy(EFields::E,GaugeFixedVariables::E);
        
        

    }
    
    
    void SetPureGauge(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "#SETTING PURE GAUGE FIELD INITIAL CONDITIONS" << std::endl;
        
        SetPureGauge(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1);
        
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
        
        
    }
    
    
    
}
