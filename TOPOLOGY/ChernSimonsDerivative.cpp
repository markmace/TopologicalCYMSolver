namespace ChernSimonsNumber{
    
    ///////////
    // NAIVE //
    ///////////
    
    // DELTA NCS ALONG THE ORIGINAL TRAJECTORY //
    DOUBLE DeltaNCsRealTime=DOUBLE(0.0);
    
    /////////////
    // COOLING //
    /////////////
    
    // DELTA NCS REAL-TIME ALONG STANDARD COOL PATH  //
    DOUBLE DeltaNCsCoolRealTime=DOUBLE(0.0);
    
    // DELTA NCS ALONG COOLING PATH //
    DOUBLE DeltaNCsCooling=DOUBLE(0.0);
    
    /////////////////
    // CALIBRATION //
    /////////////////
    
    // DELTA NCS BETWEEN STANDARD COOL PATH AND VACUUM //
    DOUBLE DeltaNCsCalibration=DOUBLE(0.0); DOUBLE DeltaNCsPreviousCalibration=DOUBLE(0.0);
    
    // DELTA NCS IN REAL TIME ALONG THE STANDARD COOLING PATH AT TIME OF CALIBRATION //
    DOUBLE DeltaNCsRealTimeCalibration=DOUBLE(0.0); DOUBLE DeltaNCsRealTimePreviousCalibration=DOUBLE(0.0);
            
    
    DOUBLE NCsDot(GaugeLinks *U,ElectricFields *E){
        
        DOUBLE EDotB;
        
        ////////////////////////////////////////
        // COMPUTE VOLUME INTEGRAL OF tr[E.B] //
        ////////////////////////////////////////
        
        // OPTION TO USE UNIMPROVED OPERATORS //
        //EDotB=UnimprovedOperators::ComputeEDotB(U,E);
        // END OPTION //
        
        // OPTION TO USE IMPROVED OPERATORS //
        EDotB=ImprovedOperators::ComputeEDotB(U,E);
        // END OPTION //
        
        
        return DOUBLE(0.125)*EDotB/SQR(PI);
    }
    
    
}