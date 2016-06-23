namespace ChernSimonsNumber{
    
    // K_0=g^2/(32 Pi^2) *LeviCitita(0,nu,rho,sigma)(F_{nu,rho}^a A_sigma^b-g/3*f^{a,b,c}*A^a_nu*A^b_rho*A^c_sigma) //
    // ONLY TAKS VACUUM PART //
    namespace VacuumEstimator{
        
        DOUBLE NCS(INT xLow, INT xHigh, INT yLow, INT yHigh, INT zLow, INT zHigh,GaugeLinks *U){
            
            // ELEMENTARY GAUGE FIELDS //
            SU_Nc_ALGEBRA_FORMAT AxM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AyM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AzM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzP[SUNcAlgebra::VectorSize];
            
            // SMEARED GAUGE FIELD //
            SU_Nc_ALGEBRA_FORMAT Ax[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Ay[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Az[SUNcAlgebra::VectorSize];
            
            // CHERN SIMONS NUMBER ESTIMATE //
            DOUBLE ChernSimonsNumber = DOUBLE(0.0);
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //CALCUALTE GAUGE FIELDS ON INCOMING AND OUTGOING LINKS //
                        
                        // OPTION TO USE NORMAL LATTICE DEFINTION //
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                        
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);
                        // END OPTION //
                        
                        // OPTION TO USE MATRIX LOGARITHM //
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                        
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                        //SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);
                        // END OPTION //
                        
                        
                        // SMEAR GAUGE LINKS AT LOCAL POINT //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            Ax[a]=DOUBLE(0.5)*(AxP[a]+AxM[a]);
                            Ay[a]=DOUBLE(0.5)*(AyP[a]+AyM[a]);
                            Az[a]=DOUBLE(0.5)*(AzP[a]+AzM[a]);
                            
                        }
                        
                        // COMPUTE LOCAL WINDING //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            for(INT b=0;b<SUNcAlgebra::VectorSize;b++){
                                for(INT c=0;c<SUNcAlgebra::VectorSize;c++){
                                    
                                    if(SUNcAlgebra::StructureFunctions::f(a,b,c) != 0){
                                        
                                        ChernSimonsNumber += DOUBLE(1.0)/(DOUBLE(16.0)*SQR(PI))*SUNcAlgebra::StructureFunctions::f(a,b,c)*Ax[a]*Ay[b]*Az[c];
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //RETURN CHERN SIMONS NUMBER ESTIMATE //
            return ChernSimonsNumber;
        }
        
        DOUBLE NCS(GaugeLinks *U){
            return NCS(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
        }
        
    }
    
}
