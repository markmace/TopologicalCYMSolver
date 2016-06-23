#ifndef __SU_Nc__GROUP_OPERATIONS__
#define __SU_Nc__GROUP_OPERATIONS__

////////////////////////////////////////////////////////////////////////////////
//																			  //
//THE MATRIX FORMAT OF THE SU(3) MATRIX U IN THE FUNDMENTAL REPRESENTATION IS //
//                                                                            //
//                          (U_0 U_3 U_6)   (U_00 U_01 U_02)                  //
//                      U=  (U_1 U_4 U_7) = (U_10 U_11 U_12)                  //
//                          (U_2 U_5 U_8)   (U_20 U_21 U_22)                  //
//                                                                            //
//                                                                            //
//																			  //
////////////////////////////////////////////////////////////////////////////////

namespace SUNcGroup{
    
    //MATRIX SIZE
    static const int MatrixSize=9;
    
    //ADJOINT SIZE
    static const int AdjointSize=64;
    
    //UNIT MATRIX
    static const SU_Nc_FUNDAMENTAL_FORMAT UnitMatrix[MatrixSize]={DOUBLE(1.0),DOUBLE(0.0),DOUBLE(0.0),DOUBLE(0.0),DOUBLE(1.0),DOUBLE(0.0),DOUBLE(0.0),DOUBLE(0.0),DOUBLE(1.0)};
    
    //INDEX CONVENTIONS
    int MatrixIndex(int i,int j){
        return i+3*j;
    }
    int AdjIndex(int a,int b){
        return a+8*b;
    }
    
    //OPERATIONS INVOLVING SU(3) MATRICES
    namespace Operations{
        
        //////////////////////
        //DETERMINANT	  //
        //////////////////////
        
        //COMPUTE THE DETERMINANT OF A COMPLEX 3x3 MATRIX
        COMPLEX Det(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return V[MatrixIndex(0,2)]*(V[MatrixIndex(1,0)]*V[MatrixIndex(2,1)]-V[MatrixIndex(1,1)]*V[MatrixIndex(2,0)])+V[MatrixIndex(0,1)]*(V[MatrixIndex(1,2)]*V[MatrixIndex(2,0)] -V[MatrixIndex(1,0)]*V[MatrixIndex(2,2)])+V[MatrixIndex(0,0)]*(V[MatrixIndex(1,1)]*V[MatrixIndex(2,2)]-V[MatrixIndex(1,2)]*V[MatrixIndex(2,1)]);
        }
        
        /////////////////////
        //INVERSE          //
        /////////////////////
        
        void Inverse(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VDagger){
            
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    
                    VDagger[MatrixIndex(i,j)]=conj(V[MatrixIndex(j,i)]);
                    
                }
            }
            
        }
        
        /////////////////////
        //TRACE OPERATIONS //
        /////////////////////
        
        //tr(V)
        COMPLEX tr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return V[MatrixIndex(0,0)]+V[MatrixIndex(1,1)]+V[MatrixIndex(2,2)];
        }
        
        DOUBLE ReTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return real(V[MatrixIndex(0,0)])+real(V[MatrixIndex(1,1)])+real(V[MatrixIndex(2,2)]);
        }
        
        DOUBLE ImTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return imag(V[MatrixIndex(0,0)])+imag(V[MatrixIndex(1,1)])+imag(V[MatrixIndex(2,2)]);
        }
        
        
        //tr(1-V)
        COMPLEX trIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(3.0)-V[MatrixIndex(0,0)]-V[MatrixIndex(1,1)]-V[MatrixIndex(2,2)];
        }
        
        DOUBLE ReTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(3.0)-real(V[MatrixIndex(0,0)])-real(V[MatrixIndex(1,1)])-real(V[MatrixIndex(2,2)]);
        }
        
        DOUBLE ImTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(3.0)-imag(V[MatrixIndex(0,0)])-imag(V[MatrixIndex(1,1)])-imag(V[MatrixIndex(2,2)]);
        }
        
        //tr(it^{a}V)
        void trIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,COMPLEX *trItV){
            
            trItV[0]= DOUBLE(0.5)*ComplexI*(V[MatrixIndex(0,1)] + V[MatrixIndex(1,0)]);
            trItV[1]=-DOUBLE(0.5)*(V[MatrixIndex(0,1)] - V[MatrixIndex(1,0)]);
            trItV[2]= DOUBLE(0.5)*ComplexI*(V[MatrixIndex(0,0)] - V[MatrixIndex(1,1)]);
            trItV[3]= DOUBLE(0.5)*ComplexI*(V[MatrixIndex(0,2)] + V[MatrixIndex(2,0)]);
            trItV[4]=-DOUBLE(0.5)*(V[MatrixIndex(0,2)] - V[MatrixIndex(2,0)]);
            trItV[5]= DOUBLE(0.5)*ComplexI*(V[MatrixIndex(1,2)] + V[MatrixIndex(2,1)]);
            trItV[6]=-DOUBLE(0.5)*(V[MatrixIndex(1,2)] - V[MatrixIndex(2,1)]);
            trItV[7]= DOUBLE(0.5)*ComplexI*(V[MatrixIndex(0,0)] + V[MatrixIndex(1,1)] - DOUBLE(2.0)*V[MatrixIndex(2,2)])/D_SQRT3;
            
        }
        
        void ReTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,1)]) + imag(V[MatrixIndex(1,0)]));
            RetrItV[1]=-DOUBLE(0.5)*(real(V[MatrixIndex(0,1)]) - real(V[MatrixIndex(1,0)]));
            RetrItV[2]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,0)]) - imag(V[MatrixIndex(1,1)]));
            RetrItV[3]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,2)]) + imag(V[MatrixIndex(2,0)]));
            RetrItV[4]=-DOUBLE(0.5)*(real(V[MatrixIndex(0,2)]) - real(V[MatrixIndex(2,0)]));
            RetrItV[5]=-DOUBLE(0.5)*(imag(V[MatrixIndex(1,2)]) + imag(V[MatrixIndex(2,1)]));
            RetrItV[6]=-DOUBLE(0.5)*(real(V[MatrixIndex(1,2)]) - real(V[MatrixIndex(2,1)]));
            RetrItV[7]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,0)]) + imag(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*imag(V[MatrixIndex(2,2)]))/D_SQRT3;
        }
        
        void ReTrIGenU(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-c*DOUBLE(0.5)*(imag(V[MatrixIndex(0,1)]) + imag(V[MatrixIndex(1,0)]));
            RetrItV[1]=-c*DOUBLE(0.5)*(real(V[MatrixIndex(0,1)]) - real(V[MatrixIndex(1,0)]));
            RetrItV[2]=-c*DOUBLE(0.5)*(imag(V[MatrixIndex(0,0)]) - imag(V[MatrixIndex(1,1)]));
            RetrItV[3]=-c*DOUBLE(0.5)*(imag(V[MatrixIndex(0,2)]) + imag(V[MatrixIndex(2,0)]));
            RetrItV[4]=-c*DOUBLE(0.5)*(real(V[MatrixIndex(0,2)]) - real(V[MatrixIndex(2,0)]));
            RetrItV[5]=-c*DOUBLE(0.5)*(imag(V[MatrixIndex(1,2)]) + imag(V[MatrixIndex(2,1)]));
            RetrItV[6]=-c*DOUBLE(0.5)*(real(V[MatrixIndex(1,2)]) - real(V[MatrixIndex(2,1)]));
            RetrItV[7]=-c*DOUBLE(0.5)*(imag(V[MatrixIndex(0,0)]) + imag(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*imag(V[MatrixIndex(2,2)]))/D_SQRT3;
        }
        
        void ImTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *ImtrItV){
            
            ImtrItV[0]= DOUBLE(0.5)*(real(V[MatrixIndex(0,1)]) + real(V[MatrixIndex(1,0)]));
            ImtrItV[1]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,1)]) - imag(V[MatrixIndex(1,0)]));
            ImtrItV[2]= DOUBLE(0.5)*(real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(1,1)]));
            ImtrItV[3]= DOUBLE(0.5)*(real(V[MatrixIndex(0,2)]) + real(V[MatrixIndex(2,0)]));
            ImtrItV[4]=-DOUBLE(0.5)*(imag(V[MatrixIndex(0,2)]) - imag(V[MatrixIndex(2,0)]));
            ImtrItV[5]= DOUBLE(0.5)*(real(V[MatrixIndex(1,2)]) + real(V[MatrixIndex(2,1)]));
            ImtrItV[6]=-DOUBLE(0.5)*(imag(V[MatrixIndex(1,2)]) - imag(V[MatrixIndex(2,1)]));
            ImtrItV[7]= DOUBLE(0.5)*(real(V[MatrixIndex(0,0)]) + real(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*real(V[MatrixIndex(2,2)]))/D_SQRT3;
            
        }
        
        
        
        //BASIC MATRIX MULTIPLICATIONS
        void UU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2){
            
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    
                    V1V2[MatrixIndex(i,j)]=0.0;
                    
                    for(int k=0;k<3;k++){
                        V1V2[MatrixIndex(i,j)]+=V1[MatrixIndex(i,k)]*V2[MatrixIndex(k,j)];
                    }
                    
                }
            }
            
        }
        
        
        void UD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2D){
            
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    
                    V1V2D[MatrixIndex(i,j)]=0.0;
                    
                    for(int k=0;k<3;k++){
                        V1V2D[MatrixIndex(i,j)]+=V1[MatrixIndex(i,k)]*conj(V2[MatrixIndex(j,k)]);
                    }
                    
                }
            }
            
        }
        
        void DU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2){
            
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    
                    V1DV2[MatrixIndex(i,j)]=0.0;
                    
                    for(int k=0;k<3;k++){
                        V1DV2[MatrixIndex(i,j)]+=conj(V1[MatrixIndex(k,i)])*V2[MatrixIndex(k,j)];
                    }
                    
                }
            }
            
        }
        
        void DD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2D){
            
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    
                    V1DV2D[MatrixIndex(i,j)]=0.0;
                    
                    for(int k=0;k<3;k++){
                        V1DV2D[MatrixIndex(i,j)]+=conj(V1[MatrixIndex(k,i)])*conj(V2[MatrixIndex(j,k)]);
                    }
                    
                }
            }
            
        }
        
        
        //CONSTRUCTION OF ADJOINT REPRESENTATION U^{adj}_{ab}=2 tr[t^{a} U_{f} t^{b} U_{f}^{\dagger}]  
        void GetAdjoint(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ADJOINT_FORMAT *VAdj){
            
            VAdj[AdjIndex(0,0)]=real(V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(0,1)]=imag(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(1,1)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(0,2)]=real(V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(0,3)]=real(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(0,4)]=imag(V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,2)]) - V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,0)]));
            
            VAdj[AdjIndex(0,5)]=real(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(0,6)]=imag(V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,2)]) - V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(0,7)]=real(V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(1,0)]=imag(V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(1,1)]=real(V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(1,2)]=imag(V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(1,3)]=imag(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(1,4)]=real(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(1,5)]=imag(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(1,6)]=real(V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,1)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(1,7)]=imag(V[MatrixIndex(1,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(1,2)]*conj(V[MatrixIndex(0,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(2,0)]=real(V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,0)]));
            
            VAdj[AdjIndex(2,1)]=imag(V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,0)])-V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,0)]));
            
            VAdj[AdjIndex(2,2)]=DOUBLE(0.5)*real(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,1)]) - V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(2,3)]=real(V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,0)]));
            
            VAdj[AdjIndex(2,4)]=imag(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,0)]));
            
            VAdj[AdjIndex(2,5)]=real(V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,1)]) - V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(2,6)]=imag(V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(2,7)]=real(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,2)]) - V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,1)]) + DOUBLE(2.0)*V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,2)]))/(DOUBLE(2.0)*D_SQRT3);
            
            VAdj[AdjIndex(3,0)]=real(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(3,1)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,1)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,0)]));
            
            VAdj[AdjIndex(3,2)]=real(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(3,3)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(3,4)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(0,0)]*conj(V[MatrixIndex(2,2)]));
            
            VAdj[AdjIndex(3,5)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(3,6)]=imag(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(0,1)]*conj(V[MatrixIndex(2,2)]));
            
            VAdj[AdjIndex(3,7)]=real(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(4,0)]=imag(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(4,1)]=real(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,1)]));
            
            VAdj[AdjIndex(4,2)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(0,1)]*conj(V[MatrixIndex(2,1)]));
            
            VAdj[AdjIndex(4,3)]=imag(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(4,4)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(4,5)]=imag(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(4,6)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,1)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,2)]));
            
            VAdj[AdjIndex(4,7)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(0,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(5,0)]=real(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(5,1)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,1)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(2,1)]));
            
            VAdj[AdjIndex(5,2)]=real(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(5,3)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(5,4)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,2)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(2,2)]));
            
            VAdj[AdjIndex(5,5)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,1)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(5,6)]=imag(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,2)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(2,2)]));
            
            VAdj[AdjIndex(5,7)]=real(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(6,0)]=imag(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(6,1)]=real(V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(6,2)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,1)]));
            
            VAdj[AdjIndex(6,3)]=imag(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(6,4)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(6,5)]=imag(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,1)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(6,6)]=real(V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,1)]) - V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,2)]));
            
            VAdj[AdjIndex(6,7)]=imag(V[MatrixIndex(2,0)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(2,1)]*conj(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(1,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,0)]=real(V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,0)]) - DOUBLE(2.0)*V[MatrixIndex(2,1)]*conj(V[MatrixIndex(2,0)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,1)]=imag(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,1)])+ V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,1)]) + DOUBLE(2.0)*V[MatrixIndex(2,1)]*conj(V[MatrixIndex(2,0)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,2)]=real(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,0)]) - V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,0)]) - V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,0)]*conj(V[MatrixIndex(2,0)]) + DOUBLE(2.0)*V[MatrixIndex(2,1)]*conj(V[MatrixIndex(2,1)]))/(DOUBLE(2.0)*D_SQRT3);
            
            VAdj[AdjIndex(7,3)]=real(V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,0)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(2,0)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,4)]=imag(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,2)]) - DOUBLE(2.0)*V[MatrixIndex(2,0)]*conj(V[MatrixIndex(2,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,5)]=real(V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,1)]) + V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(2,1)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,6)]=imag(V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,2)]) - DOUBLE(2.0)*V[MatrixIndex(2,1)]*conj(V[MatrixIndex(2,2)]))/D_SQRT3;
            
            VAdj[AdjIndex(7,7)]=real(V[MatrixIndex(0,0)]*conj(V[MatrixIndex(0,0)]) + V[MatrixIndex(0,1)]*conj(V[MatrixIndex(0,1)]) - DOUBLE(2.0)*V[MatrixIndex(0,2)]*conj(V[MatrixIndex(0,2)]) + V[MatrixIndex(1,0)]*conj(V[MatrixIndex(1,0)]) + V[MatrixIndex(1,1)]*conj(V[MatrixIndex(1,1)]) - DOUBLE(2.0)*V[MatrixIndex(1,2)]*conj(V[MatrixIndex(1,2)]) - DOUBLE(2.0)*V[MatrixIndex(2,0)]*conj(V[MatrixIndex(2,0)]) - DOUBLE(2.0)*V[MatrixIndex(2,1)]*conj(V[MatrixIndex(2,1)]) + DOUBLE(4.0)*V[MatrixIndex(2,2)]*conj(V[MatrixIndex(2,2)]))/DOUBLE(6.0);    
            
        }
        
        ///////////////////
        //UNITARITY NORM //
        ///////////////////
        
        DOUBLE UnitarityNorm(SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            SU_Nc_FUNDAMENTAL_FORMAT C[SUNcGroup::MatrixSize];
            SUNcGroup::Operations::UD(V,V,C);
            
            DOUBLE Norm=DOUBLE(0.0);
            
            for(int alpha=0;alpha<SUNcGroup::MatrixSize;alpha++){
                Norm+=SQR_ABS(C[alpha]-SUNcGroup::UnitMatrix[alpha]);
            }
            
            return sqrt(Norm);
        }
        
        
        ///////////////////////
        //DOUBLE PROJECTION  //
        ///////////////////////
        
        void ReTrIGenIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE ReTrItItV[8][8]){
            
            ReTrItItV[0][0]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(1,1)]));
            
            ReTrItItV[0][1]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,0)]) - imag(V[MatrixIndex(1,1)]));
            
            ReTrItItV[0][2]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,1)]) + real(V[MatrixIndex(1,0)]));
            
            ReTrItItV[0][3]=DOUBLE(0.25)*(-real(V[MatrixIndex(2,1)]));
            
            ReTrItItV[0][4]=DOUBLE(0.25)*(-imag(V[MatrixIndex(2,1)]));
            
            ReTrItItV[0][5]=DOUBLE(0.25)*(-real(V[MatrixIndex(2,0)]));
            
            ReTrItItV[0][6]=DOUBLE(0.25)*(-imag(V[MatrixIndex(2,0)]));
            
            ReTrItItV[0][7]=DOUBLE(0.25)*(-(real(V[MatrixIndex(0,1)])/D_SQRT3) - real(V[MatrixIndex(1,0)])/D_SQRT3);
            
            ReTrItItV[1][0]=DOUBLE(0.25)*(-imag(V[MatrixIndex(0,0)]) + imag(V[MatrixIndex(1,1)]));
            
            ReTrItItV[1][1]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(1,1)]));
            
            ReTrItItV[1][2]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,1)]) + imag(V[MatrixIndex(1,0)]));
            
            ReTrItItV[1][3]=DOUBLE(0.25)*(imag(V[MatrixIndex(2,1)]));
            
            ReTrItItV[1][4]=DOUBLE(0.25)*(-real(V[MatrixIndex(2,1)]));
            
            ReTrItItV[1][5]=DOUBLE(0.25)*(-imag(V[MatrixIndex(2,0)]));
            
            ReTrItItV[1][6]=DOUBLE(0.25)*(real(V[MatrixIndex(2,0)]));
            
            ReTrItItV[1][7]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,1)])/D_SQRT3 - imag(V[MatrixIndex(1,0)])/D_SQRT3);
            
            ReTrItItV[2][0]=DOUBLE(0.25)*(real(V[MatrixIndex(0,1)]) - real(V[MatrixIndex(1,0)]));
            
            ReTrItItV[2][1]=DOUBLE(0.25)*(-imag(V[MatrixIndex(0,1)]) - imag(V[MatrixIndex(1,0)]));
            
            ReTrItItV[2][2]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(1,1)]));
            
            ReTrItItV[2][3]=DOUBLE(0.25)*(-real(V[MatrixIndex(2,0)]));
            
            ReTrItItV[2][4]=DOUBLE(0.25)*(-imag(V[MatrixIndex(2,0)]));
            
            ReTrItItV[2][5]=DOUBLE(0.25)*(real(V[MatrixIndex(2,1)]));
            
            ReTrItItV[2][6]=DOUBLE(0.25)*(imag(V[MatrixIndex(2,1)]));
            
            ReTrItItV[2][7]=DOUBLE(0.25)*(-(real(V[MatrixIndex(0,0)])/D_SQRT3) + real(V[MatrixIndex(1,1)])/D_SQRT3);
            
            ReTrItItV[3][0]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,2)]));
            
            ReTrItItV[3][1]=DOUBLE(0.25)*(-imag(V[MatrixIndex(1,2)]));
            
            ReTrItItV[3][2]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,2)]));
            
            ReTrItItV[3][3]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(2,2)]));
            
            ReTrItItV[3][4]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,0)]) - imag(V[MatrixIndex(2,2)]));
            
            ReTrItItV[3][5]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,0)]));
            
            ReTrItItV[3][6]=DOUBLE(0.25)*(imag(V[MatrixIndex(1,0)]));
            
            ReTrItItV[3][7]=DOUBLE(0.25)*(-(real(V[MatrixIndex(0,2)])/D_SQRT3) + (2*real(V[MatrixIndex(2,0)]))/D_SQRT3);
            
            ReTrItItV[4][0]=DOUBLE(0.25)*(imag(V[MatrixIndex(1,2)]));
            
            ReTrItItV[4][1]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,2)]));
            
            ReTrItItV[4][2]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,2)]));
            
            ReTrItItV[4][3]=DOUBLE(0.25)*(-imag(V[MatrixIndex(0,0)]) + imag(V[MatrixIndex(2,2)]));
            
            ReTrItItV[4][4]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)]) - real(V[MatrixIndex(2,2)]));
            
            ReTrItItV[4][5]=DOUBLE(0.25)*(-imag(V[MatrixIndex(1,0)]));
            
            ReTrItItV[4][6]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,0)]));
            
            ReTrItItV[4][7]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,2)])/D_SQRT3 + (2*imag(V[MatrixIndex(2,0)]))/D_SQRT3);
            
            ReTrItItV[5][0]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,2)]));
            
            ReTrItItV[5][1]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,2)]));
            
            ReTrItItV[5][2]=DOUBLE(0.25)*(real(V[MatrixIndex(1,2)]));
            
            ReTrItItV[5][3]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,1)]));
            
            ReTrItItV[5][4]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,1)]));
            
            ReTrItItV[5][5]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,1)]) - real(V[MatrixIndex(2,2)]));
            
            ReTrItItV[5][6]=DOUBLE(0.25)*(imag(V[MatrixIndex(1,1)]) - imag(V[MatrixIndex(2,2)]));
            
            ReTrItItV[5][7]=DOUBLE(0.25)*(-(real(V[MatrixIndex(1,2)])/D_SQRT3) + (2*real(V[MatrixIndex(2,1)]))/D_SQRT3);
            
            ReTrItItV[6][0]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,2)]));
            
            ReTrItItV[6][1]=DOUBLE(0.25)*(real(V[MatrixIndex(0,2)]));
            
            ReTrItItV[6][2]=DOUBLE(0.25)*(-imag(V[MatrixIndex(1,2)]));
            
            ReTrItItV[6][3]=DOUBLE(0.25)*(-imag(V[MatrixIndex(0,1)]));
            
            ReTrItItV[6][4]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,1)]));
            
            ReTrItItV[6][5]=DOUBLE(0.25)*(-imag(V[MatrixIndex(1,1)]) + imag(V[MatrixIndex(2,2)]));
            
            ReTrItItV[6][6]=DOUBLE(0.25)*(-real(V[MatrixIndex(1,1)]) - real(V[MatrixIndex(2,2)]));
            
            ReTrItItV[6][7]=DOUBLE(0.25)*(imag(V[MatrixIndex(1,2)])/D_SQRT3 + (2*imag(V[MatrixIndex(2,1)]))/D_SQRT3);
            
            ReTrItItV[7][0]=DOUBLE(0.25)*(-(real(V[MatrixIndex(0,1)])/D_SQRT3) - real(V[MatrixIndex(1,0)])/D_SQRT3);
            
            ReTrItItV[7][1]=DOUBLE(0.25)*(imag(V[MatrixIndex(0,1)])/D_SQRT3 - imag(V[MatrixIndex(1,0)])/D_SQRT3);
            
            ReTrItItV[7][2]=DOUBLE(0.25)*(-(real(V[MatrixIndex(0,0)])/D_SQRT3) + real(V[MatrixIndex(1,1)])/D_SQRT3);
            
            ReTrItItV[7][3]=DOUBLE(0.25)*((2*real(V[MatrixIndex(0,2)]))/D_SQRT3 - real(V[MatrixIndex(2,0)])/D_SQRT3);
            
            ReTrItItV[7][4]=DOUBLE(0.25)*((-2*imag(V[MatrixIndex(0,2)]))/D_SQRT3 - imag(V[MatrixIndex(2,0)])/D_SQRT3);
            
            ReTrItItV[7][5]=DOUBLE(0.25)*((2*real(V[MatrixIndex(1,2)]))/D_SQRT3 - real(V[MatrixIndex(2,1)])/D_SQRT3);
            
            ReTrItItV[7][6]=DOUBLE(0.25)*((-2*imag(V[MatrixIndex(1,2)]))/D_SQRT3 - imag(V[MatrixIndex(2,1)])/D_SQRT3);
            
            ReTrItItV[7][7]=DOUBLE(0.25)*(-real(V[MatrixIndex(0,0)])/DOUBLE(3.0) - real(V[MatrixIndex(1,1)])/DOUBLE(3.0) - (DOUBLE(4.0)*real(V[MatrixIndex(2,2)]))/DOUBLE(3.0));
            
        }
        
        
        
    }
    
    namespace Extended{
        
        #include "util/CabbiboMarinariProjection.cpp"
        
        void MaxTraceProjection(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VTilde){
            
            //SET MATRIX BUFFERS
            SET_CABBIBO_MARINARI_BUFFERS();
            
            //SET INITIAL CONDITIONS
            COPY_SUNcMatrix(G,UnitMatrix);
            COPY_SUNcMatrix(VNew,V);
            
            //ITERATIVE MAX TRACE PROJECTION
            int NSteps=0;
            
            while(NSteps<7){
                
                CABBIBO_MARINARI_STEP();
                 
                NSteps++;
            }
            
            //SET FINAL RESULT
            SUNcGroup::Operations::Inverse(G,VTilde);
            
        }
        
    }
    
    namespace IO{
        
        //MATRIX TO STRING REPRESENTATION //
        std::string MatrixToString(SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            std::stringstream sstm;
            
            sstm.precision(OUTPUT_PRECISION);
            
            for(int alpha=0;alpha<MatrixSize;alpha++){
                sstm << real(V[alpha]) << " " << imag(V[alpha]) << " ";
            }
            
            return sstm.str();
            
        }
        
        //STRING TO MATRIX REPRESENTATION //
        void StringToMatrix(std::string str,SU_Nc_FUNDAMENTAL_FORMAT *V){
            
            //CONVERT TO STRING STREAM
            std::stringstream Values(str);
            
            //REAL AND IMAGINARY PARTS
            DOUBLE ReX,ImX;
            
            //GET ALL ELEMENTS
            for(int alpha=0;alpha<MatrixSize;alpha++){
                
                Values >> ReX; Values >> ImX;
                
                V[alpha]=COMPLEX(ReX,ImX);
            }
            
        }
        
        
    }
    
}

#endif
