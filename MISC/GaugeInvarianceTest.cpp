namespace GaugeInvarianceTest{
    
    
    void Check(){
        
        std::cerr << "#CHECKING GAUGE INVARIANCE" << std::endl;
        
        /////////////////////////////////////////
        // COMPUTE GAUGE INVARIANT OBSERVABLES //
        /////////////////////////////////////////
        
        Observables::Bulk::Update(GLinks::U,EFields::E);
        
        Observables::HardScales::Update(GLinks::U);
        
        std::cerr << "#ENERGY MOMENTUM TENSOR AND HARD SCALES" << std::endl;
        std::cerr << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
        
        
        DOUBLE EDotBImproved=ImprovedOperators::ComputeEDotB(GLinks::U,EFields::E);
        DOUBLE EDotBUnimproved=UnimprovedOperators::ComputeEDotB(GLinks::U,EFields::E);
        
        std::cerr << "#CHERN SIMONS TERM" << std::endl;
        std::cerr << EDotBImproved << " " << EDotBUnimproved << std::endl;
        
        /////////////////////////////////////////
        // PERFORM RANDOM GAUGE TRANSFORMATION //
        /////////////////////////////////////////
        
        std::cerr << "#SETTING RANDOM GAUGE TRANSFORMATION" << std::endl;
        GaugeTransformation::SetRandom();
        
        std::cerr << "#PERFORMING RANDOM GAUGE TRANSFORMATION" << std::endl;
        GaugeTransformation::Operations::GaugeTransformLinks();
        GaugeTransformation::Operations::GaugeTransformElectricFields();
        
        
        /////////////////////////////////////////
        // COMPUTE GAUGE INVARIANT OBSERVABLES //
        /////////////////////////////////////////
        
        Observables::Bulk::Update(GaugeFixedVariables::U,GaugeFixedVariables::E);
        
        Observables::HardScales::Update(GaugeFixedVariables::U);
        
        std::cerr << "#ENERGY MOMENTUM TENSOR AND HARD SCALES" << std::endl;
        std::cerr << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ()  << std::endl;
        
        
        DOUBLE NewEDotBImproved=ImprovedOperators::ComputeEDotB(GaugeFixedVariables::U,GaugeFixedVariables::E);
        
        DOUBLE NewEDotBUnimproved=UnimprovedOperators::ComputeEDotB(GaugeFixedVariables::U,GaugeFixedVariables::E);
        
        std::cerr << "#CHERN SIMONS TERM" << std::endl;
        std::cerr << NewEDotBImproved << " " << NewEDotBUnimproved << std::endl;
        
        std::cerr << "#DIFFERENCE IN CHERN SIMONS TERM MEASUREMENT" << std::endl;
        std::cerr << NewEDotBImproved-EDotBImproved << " " << NewEDotBUnimproved-EDotBUnimproved << std::endl;
        
        std::cerr << "#DONE CHECKING GAUGE INVARIANCE" << std::endl;
        
        
        
    }
    
}