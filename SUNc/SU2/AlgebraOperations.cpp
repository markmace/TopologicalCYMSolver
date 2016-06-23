#ifndef __SU_Nc_ALGEBRA_OPERATIONS__
#define __SU_Nc_ALGEBRA_OPERATIONS__

namespace SUNcAlgebra{
    
    //SU(2) GENERATORS ARE t^{a}=sigma_{a}/2 WHERE a=1,..3 AND sigma_{a} ARE THE PAULI MATRICES
    static const int VectorSize=3;
    
    //COMPUTES MATRIX REPRESENTATION  it^{a} alpha^a ///
    void GetMatrixFormIAlpha(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *alphaMatrix){
        
        alphaMatrix[0]=c*DOUBLE(0.5)*alpha[0];
        alphaMatrix[1]=c*DOUBLE(0.5)*alpha[1];
        alphaMatrix[2]=c*DOUBLE(0.5)*alpha[2];
        alphaMatrix[3]=0.0;

        
    }
    
    // SU(2) STRUCTURE FUNCTION //
    namespace StructureFunctions{
        
        DOUBLE f(INT a, INT b, INT c){
            return ((a-b)*(a-c)*(c-b))/DOUBLE(2.0);
        }
        
        DOUBLE d(INT I,INT j,INT k){
            return DOUBLE(0.0);
        }
        
    }
        
    //BASIC OPERATIONS
    namespace Operations{    
        
        //COMPUTES exp(c i t^{a} alpha^{a}) BY DIAGONALIZING THE MATRIX
        void MatrixIExp(DOUBLE c,DOUBLE *alpha,SU_Nc_FUNDAMENTAL_FORMAT *ExpIAlpha){
            

            DOUBLE d=sqrt(SQR(c*alpha[0])+SQR(c*alpha[1])+SQR(c*alpha[2]));
            DOUBLE arg=DOUBLE(0.5)*d;
            
            if(d!=0){
                
                DOUBLE M=c*sin(arg)/d;
                
                ExpIAlpha[0]=M*alpha[0];
                ExpIAlpha[1]=M*alpha[1];
                ExpIAlpha[2]=M*alpha[2];
                ExpIAlpha[3]=cos(arg);
            }
            else{ExpIAlpha[0]=0.0; ExpIAlpha[1]=0.0; ExpIAlpha[2]=0.0; ExpIAlpha[3]=1.0;}
            
            
        }
        
        //COMPUTES Log_{SU2}(U) BY DIAGONALIZING THE MATRIX
        //Solves exp(i c A^a t^a)=U
        void MatrixILog(DOUBLE C,SU_Nc_FUNDAMENTAL_FORMAT *U,DOUBLE *LogU){
            
            DOUBLE d=sqrt(SQR(U[0])+SQR(U[1])+SQR(U[2]));
            
            if(d!=0){
                
                DOUBLE arg=atan2(d,U[3]);
                DOUBLE M=2.0*arg/(C*d);
                
                LogU[0]=M*U[0];
                LogU[1]=M*U[1];
                LogU[2]=M*U[2];
            }
            
            else{LogU[0]=0; LogU[1]=0; LogU[2]=0;}
        }
        
        // V.E.VDagger //
        void AdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *E,SU_Nc_ALGEBRA_FORMAT *ENew){
            
            ENew[0]=E[0]*(V[0]*V[0] - V[1]*V[1] - V[2]*V[2] + V[3]*V[3]) + 2.0*(E[2]*(V[0]*V[2] - V[1]*V[3]) + E[1]*(V[0]*V[1] + V[2]*V[3]));
            ENew[1]=E[1]*(V[1]*V[1] - V[2]*V[2] - V[0]*V[0] + V[3]*V[3]) + 2.0*(E[2]*(V[1]*V[2] + V[0]*V[3]) + E[0]*(V[0]*V[1] - V[2]*V[3]));
            ENew[2]=E[2]*(V[2]*V[2] - V[0]*V[0] - V[1]*V[1] + V[3]*V[3]) + 2.0*(E[1]*(V[1]*V[2] - V[0]*V[3]) + E[0]*(V[0]*V[2] + V[1]*V[3]));
            
        }
        
        // VDagger.E.V //
        void InverseAdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *E,SU_Nc_ALGEBRA_FORMAT *ENew){
            
            ENew[0]=E[0]*(V[0]*V[0] - V[1]*V[1] - V[2]*V[2] + V[3]*V[3]) + 2.0*(E[2]*(V[0]*V[2] + V[1]*V[3]) + E[1]*(V[0]*V[1] - V[2]*V[3]));
            ENew[1]=E[1]*(V[1]*V[1] - V[2]*V[2] - V[0]*V[0] + V[3]*V[3]) + 2.0*(E[2]*(V[1]*V[2] - V[0]*V[3]) + E[0]*(V[0]*V[1] + V[2]*V[3]));
            ENew[2]=E[2]*(V[2]*V[2] - V[0]*V[0] - V[1]*V[1] + V[3]*V[3]) + 2.0*(E[1]*(V[1]*V[2] + V[0]*V[3]) + E[0]*(V[0]*V[2] - V[1]*V[3]));
        }
        
    }
    
}

#endif
