namespace InitialConditions{
    
    void SetZero(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh){
        
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    //SET ELECTRIC FIELDS TO ZERO
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            EFields::E->Get(x,y,z,mu,a)[0]=0.0;
                        }
                    }
                    
                    //SET GAUGE LINKS TO UNITY
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        COPY_SUNcMatrix(GLinks::U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                    }
                    
                }
            }
        }
    }
    
    
    void SetZero(){
        
        ///////////////////////////////
        //   SET INITIAL CONDITIONS  //
        ///////////////////////////////
        
        std::cerr << "#SETTING ZERO FIELD INITIAL CONDITIONS" << std::endl;
        
        SetZero(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1);
        
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
