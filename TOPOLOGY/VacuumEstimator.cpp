namespace ChernSimonsNumber{
    // K_0=g^2/(32 Pi^2) *LeviCitita(0,nu,rho,sigma)(F_{nu,rho}^a A_sigma^b-g/3*f^{a,b,c}*A^a_nu*A^b_rho*A^c_sigma) //
    // ONLY TAKS VACUUM PART //
    namespace VacuumEstimator{
        
        DOUBLE NCS(INT xLow, INT xHigh, INT yLow, INT yHigh, INT zLow, INT zHigh,GaugeLinks *U){
            
            // PURE GAUGE FIELDS //
            SU_Nc_ALGEBRA_FORMAT AxM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AyM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AzM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzP[SUNcAlgebra::VectorSize];
            
            // PURE GAUGE SMEARED FIELD //
            SU_Nc_ALGEBRA_FORMAT Ax[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Ay[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Az[SUNcAlgebra::VectorSize];
            
            // CHERN SIMONS NUMBER ESTIMATE //
            DOUBLE ChernSimonsNumber = DOUBLE(0.0);
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //CALCUALTE GAUGE FIELDS ON INCOMING AND OUTGOING LINKS //
                        /*
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                        
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);
                        */
                        
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                        
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);
                        
                        
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
    
    // K_0=g^2/(16 Pi^2)*LeviCitita(0,nu,rho,sigma)(A_nu^a partial_rho A_sigma^a+g/3*f^{a,b,c}*A^a_nu*A^b_rho*A^c_sigma) //
    namespace FullEstimator{
        
        DOUBLE NCS(INT xLow, INT xHigh, INT yLow, INT yHigh, INT zLow, INT zHigh,GaugeLinks *U){
            
            // PURE GAUGE FIELDS //
            SU_Nc_ALGEBRA_FORMAT AxM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AyM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyP[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AzM[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzP[SUNcAlgebra::VectorSize];
            
            // PURE GAUGE SMEARED FIELD //
            SU_Nc_ALGEBRA_FORMAT Ax[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Ay[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT Az[SUNcAlgebra::VectorSize];
            
            // PURE GAUGE DERIVATIVE FIELD //
            SU_Nc_ALGEBRA_FORMAT AxPy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxPz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyPz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyPx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzPx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzPy[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AxMy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxMz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyMz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyMx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzMx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzMy[SUNcAlgebra::VectorSize];
            
            SU_Nc_ALGEBRA_FORMAT AxMyMx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AxMzMx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyMxMy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AyMzMy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzMxMz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT AzMyMz[SUNcAlgebra::VectorSize];

            
            SU_Nc_ALGEBRA_FORMAT DzAx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT DyAx[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT DxAy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT DzAy[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT DxAz[SUNcAlgebra::VectorSize];
            SU_Nc_ALGEBRA_FORMAT DyAz[SUNcAlgebra::VectorSize];

            
            // CHERN SIMONS NUMBER ESTIMATE //
            DOUBLE ChernSimonsNumber = DOUBLE(0.0);
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        //CALCUALTE GAUGE FIELDS ON INCOMING AND OUTGOING LINKS //
                        /*
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                         
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                         SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);
                         */
                        
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,0),AxP);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,1),AyP);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z,2),AzP);
                        
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y,z,0),AxM);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y-1,z,1),AyM);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z-1,2),AzM);

                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y+1,z,0),AxPy);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z+1,0),AxPz);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y,z+1,1),AyPz);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x+1,y,z,1),AyPx);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x+1,y,z,2),AzPx);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y+1,z,2),AzPy);
                        
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y-1,z,0),AxMyMx);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y,z-1,0),AxMzMx);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y-1,z-1,1),AyMzMy);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y-1,z,1),AyMxMy);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x-1,y,z-1,2),AzMxMz);
                        SUNcAlgebra::Operations::MatrixILog(DOUBLE(2.0),U->Get(x,y-1,z-1,2),AzMyMz);
            
                        
                        // SMEAR GAUGE LINKS AT LOCAL POINT //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            Ax[a]=DOUBLE(0.5)*(AxP[a]+AxM[a]);
                            Ay[a]=DOUBLE(0.5)*(AyP[a]+AyM[a]);
                            Az[a]=DOUBLE(0.5)*(AzP[a]+AzM[a]);
                            
                            DyAx[a]=DOUBLE(0.25)*(AxPy[a]+AxP[a]-AxMy[a]-AxMyMx[a]);
                            DzAx[a]=DOUBLE(0.25)*(AxPz[a]+AxP[a]-AxMz[a]-AxMzMx[a]);
                            DzAy[a]=DOUBLE(0.25)*(AyPz[a]+AyP[a]-AyMz[a]-AyMzMy[a]);
                            DxAy[a]=DOUBLE(0.25)*(AyPx[a]+AyP[a]-AyMx[a]-AyMxMy[a]);
                            DxAz[a]=DOUBLE(0.25)*(AzPx[a]+AzP[a]-AzMx[a]-AzMxMz[a]);
                            DyAz[a]=DOUBLE(0.25)*(AzPy[a]+AzP[a]-AzMy[a]-AzMyMz[a]);

                        }
                        
                        // COMPUTE LOCAL WINDING //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            // NON VACUUM PART
                            ChernSimonsNumber += DOUBLE(1.0)/(DOUBLE(16.0)*SQR(PI))*(Ax[a]*(DyAz[a]-DzAy[a])+Ay[a]*(DzAx[a]-DxAz[a])+Az[a]*(DxAy[a]-DyAx[a]));
                            
                            for(INT b=0;b<SUNcAlgebra::VectorSize;b++){
                                for(INT c=0;c<SUNcAlgebra::VectorSize;c++){
                                    
                                    
                                    // VACUUM PART //
                                    
                                    // IS THIS OFF BY A MINUS SIGN OR FACTOR!!???!!
                                    if(SUNcAlgebra::StructureFunctions::f(a,b,c) != 0){
                                        
                                        ChernSimonsNumber += DOUBLE(-1.0)/(DOUBLE(8.0)*SQR(PI))*SUNcAlgebra::StructureFunctions::f(a,b,c)*Ax[a]*Ay[b]*Az[c];
                                        
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
