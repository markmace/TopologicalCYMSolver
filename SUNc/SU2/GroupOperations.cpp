#ifndef __SU_Nc__GROUP_OPERATIONS__
#define __SU_Nc__GROUP_OPERATIONS__

////////////////////////////////////////////////////////////////////////////////
//																			  //
//THE MATRIX FORMAT OF THE SU(2) MATRIX U IN THE FUNDMENTAL REPRESENTATION IS //
//				U=u[3]*ID+u[i] I Sigma[i] i=0..(Nc^2-2)						  //
//																			  //
////////////////////////////////////////////////////////////////////////////////

namespace SUNcGroup{
    
    //MATRIX SIZE
    static const int MatrixSize=4;
    
    //ADJOINT SIZE
    static const int AdjointSize=9;
    
    //UNIT MATRIX
    static const SU_Nc_FUNDAMENTAL_FORMAT UnitMatrix[MatrixSize]={DOUBLE(0.0),DOUBLE(0.0),DOUBLE(0.0),DOUBLE(1.0)};
    
    //INDEX CONVENTIONS
    int AdjIndex(int a,int b){
        return a+3*b;
    }
    
    //OPERATIONS INVOLVING SU(2) MATRICES
    namespace Operations{
        
      	COMPLEX Det(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return SQR(V[0])+SQR(V[1])+SQR(V[2])+SQR(V[3]);
        }
        
        /////////////////////
        //INVERSE          //
        /////////////////////
        
        void Inverse(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VDagger){
            
            VDagger[0]=-V[0];
            VDagger[1]=-V[1];
            VDagger[2]=-V[2];
            VDagger[3]=V[3];
            
        }
        
        /////////////////////
        //TRACE OPERATIONS //
        /////////////////////
        
        //tr(V)
        COMPLEX tr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(2.0)*V[3];
        }
        
        DOUBLE ReTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(2.0)*V[3];
        }
        
        DOUBLE ImTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(0.0);
        }
        
        
        //tr(1-V)
        COMPLEX trIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(2.0)*(DOUBLE(1.0)-V[3]);
        }
        
        DOUBLE ReTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(2.0)*(DOUBLE(1.0)-V[3]);
        }
        
        DOUBLE ImTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(0.0);
        }
        
        void ReTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-V[0];
            RetrItV[1]=-V[1];
            RetrItV[2]=-V[2];
            
        }
        
        void ReTrIGenU(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-c*V[0];
            RetrItV[1]=-c*V[1];
            RetrItV[2]=-c*V[2];
            
        }
        
        void ImTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *ImtrItV){
            
            ImtrItV[0]=0.0;
            ImtrItV[0]=0.0;
            ImtrItV[0]=0.0;
            
        }
        
        //BASIC MATRIX MULTIPLICATIONS
        void UU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2){
            
            V1V2[0]=V1[3]*V2[0]+V1[0]*V2[3]-(V1[1]*V2[2]-V1[2]*V2[1]);
            V1V2[1]=V1[3]*V2[1]+V1[1]*V2[3]+(V1[0]*V2[2]-V1[2]*V2[0]);
            V1V2[2]=V1[3]*V2[2]+V1[2]*V2[3]-(V1[0]*V2[1]-V1[1]*V2[0]);
            V1V2[3]=V1[3]*V2[3]-V1[0]*V2[0]-V1[1]*V2[1]-V1[2]*V2[2];
            
        }
        
        
        void UD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2D){
            
            V1V2D[0]=-V1[3]*V2[0]+V1[0]*V2[3]+(V1[1]*V2[2]-V1[2]*V2[1]);
            V1V2D[1]=-V1[3]*V2[1]+V1[1]*V2[3]-(V1[0]*V2[2]-V1[2]*V2[0]);
            V1V2D[2]=-V1[3]*V2[2]+V1[2]*V2[3]+(V1[0]*V2[1]-V1[1]*V2[0]);
            V1V2D[3]=V1[3]*V2[3]+V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
            
        }
        
        void DU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2){
            
            V1DV2[0]=V1[3]*V2[0]-V1[0]*V2[3]+(V1[1]*V2[2]-V1[2]*V2[1]);
            V1DV2[1]=V1[3]*V2[1]-V1[1]*V2[3]-(V1[0]*V2[2]-V1[2]*V2[0]);
            V1DV2[2]=V1[3]*V2[2]-V1[2]*V2[3]+(V1[0]*V2[1]-V1[1]*V2[0]);
            V1DV2[3]=V1[3]*V2[3]+V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
            
        }
        
        void DD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2D){
            
            V1DV2D[0]=-V1[3]*V2[0]-V1[0]*V2[3]-(V1[1]*V2[2]-V1[2]*V2[1]);
            V1DV2D[1]=-V1[3]*V2[1]-V1[1]*V2[3]+(V1[0]*V2[2]-V1[2]*V2[0]);
            V1DV2D[2]=-V1[3]*V2[2]-V1[2]*V2[3]-(V1[0]*V2[1]-V1[1]*V2[0]);
            V1DV2D[3]=V1[3]*V2[3]-V1[0]*V2[0]-V1[1]*V2[1]-V1[2]*V2[2];
            
        }
        
        
        //CONSTRUCTION OF ADJOINT REPRESENTATION U^{adj}_{ab}=2 tr[t^{a} U_{f} t^{b} U_{f}^{\dagger}]  
        void GetAdjoint(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ADJOINT_FORMAT *VAdj){
            
            VAdj[AdjIndex(0,0)]= SQR(V[0])+ SQR(V[3]) - SQR(V[1]) - SQR(V[2]) ;
            VAdj[AdjIndex(0,1)]= DOUBLE(2.0)*(V[0]*V[1]+V[2]*V[3]);
            VAdj[AdjIndex(0,2)]= DOUBLE(2.0)*(V[0]*V[2]-V[1]*V[3]);
            VAdj[AdjIndex(1,0)]= DOUBLE(2.0)*(V[0]*V[1]-V[2]*V[3]);
            VAdj[AdjIndex(1,1)]= SQR(V[1])+ SQR(V[3]) - SQR(V[0]) - SQR(V[2]) ;
            VAdj[AdjIndex(1,2)]= DOUBLE(2.0)*(V[1]*V[2]+V[0]*V[3]);
            VAdj[AdjIndex(2,0)]= DOUBLE(2.0)*(V[0]*V[2]+V[1]*V[3]);
            VAdj[AdjIndex(2,1)]= DOUBLE(2.0)*(V[1]*V[2]-V[0]*V[3]);
            VAdj[AdjIndex(2,2)]= SQR(V[2])+ SQR(V[3]) - SQR(V[0]) - SQR(V[1]) ;
            
        }
        
        /////////////////////////
        //PROJECTION to SPHERE //
        /////////////////////////
        
        void SphereProj(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *SphereProjV){
            
            SphereProjV[0]=-V[0];
            SphereProjV[1]=-V[1];
            SphereProjV[2]=-V[2];
            SphereProjV[3]=V[3];
        }
        
        ///////////////////
        //UNITARITY NORM //
        ///////////////////
        
        DOUBLE UnitarityNorm(SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            SU_Nc_FUNDAMENTAL_FORMAT C[SUNcGroup::MatrixSize];
            SUNcGroup::Operations::UD(V,V,C);
            
            DOUBLE Norm=DOUBLE(0.0);
            
            for(int alpha=0;alpha<SUNcGroup::MatrixSize;alpha++){
                Norm+=SQR(C[alpha]-SUNcGroup::UnitMatrix[alpha]);
            }
            
            return sqrt(Norm);
        }
        
        ///////////////////////
        //DOUBLE PROJECTION  //
        ///////////////////////
        
        void ReTrIGenIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE ReTrItItV[3][3]){
            
            ReTrItItV[0][0]=-0.5*V[3];
            
            ReTrItItV[0][1]=0.5*V[2];
            
            ReTrItItV[0][2]=-0.5*V[1];
            
            ReTrItItV[1][0]=-0.5*V[2];
            
            ReTrItItV[1][1]=-0.5*V[3];
            
            ReTrItItV[1][2]=0.5*V[0];
            
            ReTrItItV[2][0]=0.5*V[1];
            
            ReTrItItV[2][1]=-0.5*V[0];
            
            ReTrItItV[2][2]=-0.5*V[3];

            
        }

        
    }
    
    namespace Extended{
        
        // COMPUTE UNITARIZATION UTilde=U/Sqrt(U^{dagger} U) //
        void MaxTraceProjection(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VTilde){
            
            DOUBLE Norm=sqrt(SQR(V[0])+SQR(V[1])+SQR(V[2])+SQR(V[3]));
            
            VTilde[0]=V[0]/Norm;
            VTilde[1]=V[1]/Norm;
            VTilde[2]=V[2]/Norm;
            VTilde[3]=V[3]/Norm;
            
        }
        
        // COMPUTE UNITARIZATION UTilde=U/Sqrt(U^{dagger} U) -- IN-PLACE OPTION //
        void MaxTraceProjection(SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            DOUBLE Norm=sqrt(SQR(V[0])+SQR(V[1])+SQR(V[2])+SQR(V[3]));
            
            V[0]=V[0]/Norm;
            V[1]=V[1]/Norm;
            V[2]=V[2]/Norm;
            V[3]=V[3]/Norm;
            
        }
    }
    
    namespace IO{
        
        //MATRIX TO STRING REPRESENTATION //
        std::string MatrixToString(SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            std::stringstream sstm;
            
            sstm.precision(OUTPUT_PRECISION);
            
            for(int alpha=0;alpha<MatrixSize;alpha++){
                sstm << V[alpha] << " ";
            }
            
            return sstm.str();
            
        }
        
        //STRING TO MATRIX REPRESENTATION //
        void StringToMatrix(std::string str,SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            //CONVERT TO STRING STREAM
            std::stringstream Values(str);
            
            //REAL AND IMAGINARY PARTS
            DOUBLE ReX;
            
            //GET ALL ELEMENTS
            for(int alpha=0;alpha<MatrixSize;alpha++){
                
                Values >> ReX;
                
                V[alpha]=DOUBLE(ReX);
            }
            
        }
    }
    
}

#endif
