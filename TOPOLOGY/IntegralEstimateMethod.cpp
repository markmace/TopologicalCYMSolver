namespace ChernSimonsNumber{
    
    namespace IntegralEstimateMethod{
        
        DOUBLE NCS(INT xLow, INT xHigh, INT yLow, INT yHigh, INT zLow, INT zHigh,GaugeTransformations *G){
            
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
            
            // SU(Nc) MATRIX BUFFER //
            SU_Nc_FUNDAMENTAL_FORMAT MatrixBuffer[SUNcGroup::MatrixSize];
            
            // CHERN SIMONS NUMBER ESTIMATE //
            DOUBLE ChernSimonsNumber = DOUBLE(0.0);
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        //DEFINED AS 1/2*ReTr(S(x+i)^\dagger (-i \lambda^a) S(x)), we compute 1/2*ReTr( (-i \lambda^a) S(x) S(x+i)^\dagger)//
                        //CALCUALTE GAUGE FIELDS ON INCOMING AND OUTGOING LINKS //
                        SUNcGroup::Operations::UD(G->Get(x,y,z),G->Get(x+1,y,z),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AxP);
                        
                        SUNcGroup::Operations::UD(G->Get(x-1,y,z),G->Get(x,y,z),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AxM);
                        
                        SUNcGroup::Operations::UD(G->Get(x,y,z),G->Get(x,y+1,z),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AyP);
                        
                        SUNcGroup::Operations::UD(G->Get(x,y-1,z),G->Get(x,y,z),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AyM);
                        
                        SUNcGroup::Operations::UD(G->Get(x,y,z),G->Get(x,y,z+1),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AzP);
                        
                        SUNcGroup::Operations::UD(G->Get(x,y,z-1),G->Get(x,y,z),MatrixBuffer);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),MatrixBuffer,AzM);
                        
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
        
        DOUBLE NCS(GaugeLinks *U,GaugeTransformations *G){
            return NCS(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,G);
        }
        
        DOUBLE NCS(GaugeTransformations *G){
            return NCS(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,G);
        }
        
        DOUBLE NCS(){
            return NCS(GaugeTransformation::G);
        }
        
    }

    
}


