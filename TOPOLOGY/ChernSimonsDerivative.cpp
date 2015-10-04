namespace ChernSimonsNumber{
    
    DOUBLE DeltaNCsRealTime=DOUBLE(0.0);  DOUBLE DeltaNCsCoolRealTime=DOUBLE(0.0);  DOUBLE DeltaNCsCooling=DOUBLE(0.0); DOUBLE DeltaNCsPreviousCooling=DOUBLE(0.0); DOUBLE DeltaNCsVacuum=DOUBLE(0.0);
            
    // COMPUTE DERIVATIVE OF THE CHERN SIMONS NUMBER //
    DOUBLE NCsDot(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        /*
        // RESET //
        DOUBLE Result=DOUBLE(0.0);
        
        DOUBLE CumESqr=DOUBLE(0.0);
        DOUBLE CumBSqr=DOUBLE(0.0);
        
        #pragma omp parallel
        {
            
            //BUFFERS FOR AVERAGE FIELD STRENGTH
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            // SET CONSTANTS //
            DOUBLE cS[Lattice::Dimension];
            
            for(INT mu=0;mu<Lattice::Dimension;mu++){
                cS[mu]=Lattice::aScale*SQR(U->a[mu])/U->aCube;
            }
        
            //UPDATE AT ALL SITES //
            #pragma omp for reduction( + : Result,CumESqr,CumBSqr)
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        // COMPUTE AVERAGE E and B FIELDS //
                        COMPUTE_AVG_FIELD_STRENGTH(x,y,z);
                        
                        // UPDATE E.B //
                        Result+=cS[0]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B0Loc);
                        Result+=cS[1]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B1Loc);
                        Result+=cS[2]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B2Loc);

                        
                        CumBSqr+=B0SqrLoc+B1SqrLoc+B2SqrLoc;
                        CumESqr+=E0SqrLoc+E1SqrLoc+E2SqrLoc;
                        

                    }
                }
            }
        
            
        }// END PARALLEL

        
        // RETURN NORMALIZED RESULT //
        return DOUBLE(0.125)*Result/SQR(PI);
        */
        
        
        DOUBLE EdB=ImprovedOperators::Compute(U,E);
        return DOUBLE(0.125)*EdB/SQR(PI);


        
    }
    
    DOUBLE NCsDot(GaugeLinks *U,ElectricFields *E){
        
        return NCsDot(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
    }
    
    
    
}