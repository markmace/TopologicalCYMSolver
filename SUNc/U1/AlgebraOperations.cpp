#ifndef __SU_Nc_ALGEBRA_OPERATIONS__
#define __SU_Nc_ALGEBRA_OPERATIONS__

namespace SUNcAlgebra{
    
    //U(1) GENERATORS ARE t=1/sqrt(2)
    static const int VectorSize=1;
    
    //COMPUTES MATRIX REPRESENTATION  it^{a}alpha^a ///
    void GetMatrixFormIAlpha(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *alphaMatrix){
        
        alphaMatrix[0]=c*COMPLEX(0.0,0.5)*DOUBLE(M_SQRT2)*alpha[0];
        
    }
        
    //BASIC OPERATIONS
    namespace Operations{    
        
        //COMPUTES exp( c i t^{a} alpha^{a}) BY DIAGONALIZING THE MATRIX
        void MatrixIExp(DOUBLE c,DOUBLE *alpha,SU_Nc_FUNDAMENTAL_FORMAT *ExpIAlpha){
            
            ExpIAlpha[0]=COMPLEX (cos(DOUBLE(0.5)*c*M_SQRT2*alpha[0]),sin(DOUBLE(0.5)*c*M_SQRT2*alpha[0]));
            
        }
        
        // COMPUTE gamma=-i[alpha,beta] //
        void MinusILieBracket(SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_ALGEBRA_FORMAT *beta,SU_Nc_ALGEBRA_FORMAT *gamma){
            gamma[0]=DOUBLE(0.0);
        }
        
    }
    
}

#endif
