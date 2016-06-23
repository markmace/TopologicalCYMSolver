#ifndef __SU_NC_ALGEBRA_OPERATIONS__
#define __SU_NC_ALGEBRA_OPERATIONS__

namespace SUNcAlgebra{
    
    //SU(3) GENERATORS ARE t^{a}=lambda_{a}/2 WHERE lambda_{a} (a=1,..,8) ARE THE GELL-MANN MATRICES
    static const int VectorSize=8;
    
    //COMPUTES MATRIX REPRESENTATION  it^{a} alpha^a ///
    void GetMatrixFormIAlpha(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *alphaMatrix){
        
        alphaMatrix[SUNcGroup::MatrixIndex(0,0)]= c*COMPLEX(0.0,0.5)*(alpha[2]+alpha[7]/D_SQRT3);
        alphaMatrix[SUNcGroup::MatrixIndex(1,1)]=-c*COMPLEX(0.0,0.5)*(alpha[2]-alpha[7]/D_SQRT3);
        alphaMatrix[SUNcGroup::MatrixIndex(2,2)]=-c*COMPLEX(0.0,alpha[7]/D_SQRT3);
        
        alphaMatrix[SUNcGroup::MatrixIndex(1,0)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[0],alpha[1]);
        alphaMatrix[SUNcGroup::MatrixIndex(2,0)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[3],alpha[4]);
        alphaMatrix[SUNcGroup::MatrixIndex(2,1)]=c*COMPLEX(0.0,0.5)*COMPLEX(alpha[5],alpha[6]);
        
        alphaMatrix[SUNcGroup::MatrixIndex(0,1)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(1,0)]);
        alphaMatrix[SUNcGroup::MatrixIndex(0,2)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(2,0)]);
        alphaMatrix[SUNcGroup::MatrixIndex(1,2)]=conj(alphaMatrix[SUNcGroup::MatrixIndex(2,1)]);
        
        
    }
    
    // SU(3) STRUCTURE FUNCTION //
    namespace StructureFunctions{
        // ANTI-SYMMETRIC TENSOR
        DOUBLE f(INT a, INT b, INT c){
            
            if(a==b || b==c || a==c){
                return 0.0;
            }
            // f^{123}=1
            else if((a==1 || b==1 || c==1) && (a==2 || b==2 || c==2) && (a==3 || b==3 || c==3)){
                return ((a-b)*(a-c)*(c-b))/DOUBLE(2.0);
            }
            // f^{147}=1/2
            else if((a==1 || b==1 || c==1) && (a==4 || b==4 || c==4) && (a==7 || b==7 || c==7)){
                return DOUBLE(0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(54.0);
            }
            // f^{246}=1/2
            else if((a==2 || b==2 || c==2) && (a==4 || b==4 || c==4) && (a==6 || b==6 || c==6)){
                return DOUBLE(0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(16.0);
            }
            // f^{257}=1/2
            else if((a==2 || b==2 || c==2) && (a==5 || b==5 || c==5) && (a==7 || b==7 || c==7)){
                return DOUBLE(0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(30.0);
            }
            // f^{345}=1/2
            else if((a==3 || b==3 || c==3) && (a==4 || b==4 || c==4) && (a==5 || b==5 || c==5)){
                return DOUBLE(0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(2.0);
            }
            // f^{367}=-1/2
            else if((a==3 || b==3 || c==3) && (a==6 || b==6 || c==6) && (a==7 || b==7 || c==7)){
                return DOUBLE(-0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(12.0);
            }
            // f^{156}=-1/2
            else if((a==1 || b==1 || c==1) && (a==5 || b==5 || c==5) && (a==6 || b==6 || c==6)){
                return DOUBLE(-0.5)*((a-b)*(a-c)*(c-b))/DOUBLE(20.0);
            }
            // f^{458}=sqrt(3)/2
            else if((a==4 || b==4 || c==4) && (a==5 || b==5 || c==5) && (a==8 || b==8 || c==8)){
                return DOUBLE(0.5*D_SQRT3)*((a-b)*(a-c)*(c-b))/DOUBLE(12.0);
            }
            // f^{678}=sqrt(3)/2
            else if((a==6 || b==6 || c==6) && (a==7 || b==7 || c==7) && (a==8 || b==8 || c==8)){
                return DOUBLE(0.5*D_SQRT3)*((a-b)*(a-c)*(c-b))/DOUBLE(2.0);
            }
            else{
                return 0.0;
            }
        }
        
        /*
        // NEEDS TO BE ADDED //
        DOUBLE d(INT I,INT j,INT k){
            return DOUBLE(0.0);
        }
         */
        
    }
    
    //BASIC OPERATIONS
    namespace Operations{    
        
        //COMPUTES MATRIX REPRESENTATION  t^{a} alpha^a ///
        void MatrixForm(SU_Nc_ALGEBRA_FORMAT *alpha,COMPLEX A[3][3]){
            
            A[0][0]= DOUBLE(0.5)*(alpha[2]+alpha[7]/D_SQRT3);
            A[1][1]=-DOUBLE(0.5)*(alpha[2]-alpha[7]/D_SQRT3);
            A[2][2]=-alpha[7]/D_SQRT3;
            
            A[1][0]=DOUBLE(0.5)*COMPLEX(alpha[0],alpha[1]);
            A[2][0]=DOUBLE(0.5)*COMPLEX(alpha[3],alpha[4]);
            A[2][1]=DOUBLE(0.5)*COMPLEX(alpha[5],alpha[6]);
            
            A[0][1]=conj(A[1][0]);
            A[0][2]=conj(A[2][0]);
            A[1][2]=conj(A[2][1]);
            
        }

        
        //CHANGES SUNC FUNDAMENTAL MATRIX FROM U[i=0,..,8] to U[i=0,..,2][j=0,..,2] ///
        void ChangeFundamentalForm1to2(SU_Nc_FUNDAMENTAL_FORMAT *U1,SU_Nc_FUNDAMENTAL_FORMAT U2[3][3]){
            
            U2[0][0]=U1[SUNcGroup::MatrixIndex(0,0)];
            U2[0][1]=U1[SUNcGroup::MatrixIndex(0,1)];
            U2[0][2]=U1[SUNcGroup::MatrixIndex(0,2)];
            
            U2[1][0]=U1[SUNcGroup::MatrixIndex(1,0)];
            U2[1][1]=U1[SUNcGroup::MatrixIndex(1,1)];
            U2[1][2]=U1[SUNcGroup::MatrixIndex(1,2)];
            
            U2[2][0]=U1[SUNcGroup::MatrixIndex(2,0)];
            U2[2][1]=U1[SUNcGroup::MatrixIndex(2,1)];
            U2[2][2]=U1[SUNcGroup::MatrixIndex(2,2)];
            
        }
        
        //CHANGES SUNC FUNDAMENTAL MATRIX FROM U[i=0,..,2][j=0,..,2] to U[i=0,..,8] ///
        void ChangeFundamentalForm2to1(SU_Nc_FUNDAMENTAL_FORMAT U2[3][3],SU_Nc_FUNDAMENTAL_FORMAT U1[SUNcGroup::MatrixSize]){
            
            U1[SUNcGroup::MatrixIndex(0,0)]=U2[0][0];
            U1[SUNcGroup::MatrixIndex(0,1)]=U2[0][1];
            U1[SUNcGroup::MatrixIndex(0,2)]=U2[0][2];
            
            U1[SUNcGroup::MatrixIndex(1,0)]=U2[1][0];
            U1[SUNcGroup::MatrixIndex(1,1)]=U2[1][1];
            U1[SUNcGroup::MatrixIndex(1,2)]=U2[1][2];
            
            U1[SUNcGroup::MatrixIndex(2,0)]=U2[2][0];
            U1[SUNcGroup::MatrixIndex(2,1)]=U2[2][1];
            U1[SUNcGroup::MatrixIndex(2,2)]=U2[2][2];
            
        }
        
        //COMPUTES exp( c i t^{a} alpha^{a}) BY DIAGONALIZING THE MATRIX
        void MatrixIExp(DOUBLE c,SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_FUNDAMENTAL_FORMAT *ExpIAlpha){
            
            //GET MATRIX REPRESENTATION
            COMPLEX A[3][3];
            MatrixForm(alpha,A);
            
            //COMPUTE EIGENSYSTEM
            COMPLEX U[3][3];
            DOUBLE eV[3];
            
            int SUCCESS=Diagonalization3x3::Eigensystem(A,U,eV);
            
            if(SUCCESS==-1){
                std::cerr << "DIAGONALIZATION FAILED!" << std::endl;
                exit(0);
            }
            
            //COMPUTE MATRIX EXPONENTIAL
            
            //SET TO ZERO
            for(int i=0;i<3;i++){
                for(int k=0;k<3;k++){
                    ExpIAlpha[SUNcGroup::MatrixIndex(i,k)]=DOUBLE(0.0);
                }
            }
            
            //PERFORM EXPONENTIATION
            for(int j=0;j<3;j++){
                
                //COMPUTE EXPONENTIAL OF I TIMES THE EIGENVALUE TIMES THE CONSTANT FACTOR c
                COMPLEX ExpIEv=COMPLEX(cos(c*eV[j]),sin(c*eV[j]));
                
                //CONTRACT WITH MATRIX
                for(int i=0;i<3;i++){
                    for(int k=0;k<3;k++){
                        ExpIAlpha[SUNcGroup::MatrixIndex(i,k)]+=U[i][j]*ExpIEv*conj(U[k][j]);
                    }
                }
                
            }
            
        }
        
        //COMPUTES Log_{SU3}(U) BY DIAGONALIZING THE MATRIX
        //Solves exp(i c A^a t^a)=U
        void MatrixILog(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *U,SU_Nc_ALGEBRA_FORMAT *MinusILogU){
            
            // AS TO NOT OVERWRITE ORIGINAL MATRIX
            SU_Nc_FUNDAMENTAL_FORMAT BUFFER[SUNcGroup::MatrixSize];
            COPY_SUNcMatrix(BUFFER,U);
            
            //COMPUTE EIGENSYSTEM
            COMPLEX eVals[3];
            COMPLEX eVecs[9];
            
            // WILL WRITE OVER BUFFER
            int SUCCESS=ComplexEigensystem::GetComplexEigensystem(BUFFER,eVals,eVecs);
            
            if(SUCCESS==-1){
                std::cerr << "DIAGONALIZATION FAILED!" << std::endl;
                exit(0);
            }
            
            //INITIALLY SET TO ZERO
            for(int i=0;i<SUNcAlgebra::VectorSize;i++){
                MinusILogU[i]=DOUBLE(0.0);
            }
            
            SU_Nc_FUNDAMENTAL_FORMAT LogUMat[SUNcGroup::MatrixSize];
         
            //SET DIAGONAL LOG MATRIX //
            //PERFORM EXPONENTIATION
            for(INT j=0;j<3;j++){
                
                //COMPUTE INVERSE OF EXPONETIAL FORM OF EIGENVALUE TIMES THE CONSTANT FACTOR c
                // LOG(z)=i ARG(z)
                COMPLEX IArgEv=c*arg(eVals[j])*COMPLEX(0.0,1.0);
                
                //CONTRACT WITH MATRIX
                for(INT i=0;i<3;i++){
                    for(INT k=0;k<3;k++){
                        LogUMat[SUNcGroup::MatrixIndex(i,k)]+=eVecs[SUNcGroup::MatrixIndex(i,j)]*IArgEv*conj(eVecs[SUNcGroup::MatrixIndex(k,j)]);
                        
                    }
                }
                
            }
            
            // -I*LOG(U)
            SUNcGroup::Operations::ReTrIGenU(-2.0,LogUMat,MinusILogU);
            
            
        }
        
        // COMPUTE gamma=-i[alpha,beta] //
        void MinusILieBracket(SU_Nc_ALGEBRA_FORMAT *alpha,SU_Nc_ALGEBRA_FORMAT *beta,SU_Nc_ALGEBRA_FORMAT *gamma){
            
            gamma[0]=-alpha[2]*beta[1]+alpha[1]*beta[2]-DOUBLE(0.5)*alpha[6]*beta[3]+DOUBLE(0.5)*alpha[5]*beta[4]-DOUBLE(0.5)*alpha[4]*beta[5]+DOUBLE(0.5)*alpha[3]*beta[6];
            
            gamma[1]=alpha[2]*beta[0]-alpha[0]*beta[2]-DOUBLE(0.5)*alpha[5]*beta[3]-DOUBLE(0.5)*alpha[6]*beta[4]+DOUBLE(0.5)*alpha[3]*beta[5]+DOUBLE(0.5)*alpha[4]*beta[6];
            
            gamma[2]=-alpha[1]*beta[0]+alpha[0]*beta[1]-DOUBLE(0.5)*alpha[4]*beta[3]+DOUBLE(0.5)*alpha[3]*beta[4]+DOUBLE(0.5)*alpha[6]*beta[5]-DOUBLE(0.5)*alpha[5]*beta[6];
            
            gamma[3]= DOUBLE(0.5)*alpha[6]*beta[0]+DOUBLE(0.5)*alpha[5]*beta[1]+DOUBLE(0.5)*alpha[4]*beta[2]-DOUBLE(0.5)*alpha[2]*beta[4]-DOUBLE(0.5)*D_SQRT3*alpha[7]*beta[4]-DOUBLE(0.5)*alpha[1]*beta[5]-DOUBLE(0.5)*alpha[0]*beta[6]+DOUBLE(0.5)*D_SQRT3*alpha[4]*beta[7];
            
            gamma[4]=-DOUBLE(0.5)*alpha[5]*beta[0]+DOUBLE(0.5)*alpha[6]*beta[1]-DOUBLE(0.5)*alpha[3]*beta[2]+DOUBLE(0.5)*alpha[2]*beta[3]+DOUBLE(0.5)*D_SQRT3*alpha[7]*beta[3]+DOUBLE(0.5)*alpha[0]*beta[5]-DOUBLE(0.5)*alpha[1]*beta[6]-DOUBLE(0.5)*D_SQRT3*alpha[3]*beta[7];
            
            gamma[5]=DOUBLE(0.5)*alpha[4]*beta[0]-DOUBLE(0.5)*alpha[3]*beta[1]-DOUBLE(0.5)*alpha[6]*beta[2]+DOUBLE(0.5)*alpha[1]*beta[3]-DOUBLE(0.5)*alpha[0]*beta[4]+DOUBLE(0.5)*alpha[2]*beta[6]-DOUBLE(0.5)*D_SQRT3*alpha[7]*beta[6]+DOUBLE(0.5)*D_SQRT3*alpha[6]*beta[7];
            
            gamma[6]=-DOUBLE(0.5)*alpha[3]*beta[0]-DOUBLE(0.5)*alpha[4]*beta[1]+DOUBLE(0.5)*alpha[5]*beta[2]+DOUBLE(0.5)*alpha[0]*beta[3]+DOUBLE(0.5)*alpha[1]*beta[4]-DOUBLE(0.5)*alpha[2]*beta[5]+DOUBLE(0.5)*D_SQRT3*alpha[7]*beta[5]-DOUBLE(0.5)*D_SQRT3*alpha[5]*beta[7];
            
            gamma[7]=-DOUBLE(0.5)*D_SQRT3*alpha[4]*beta[3]+DOUBLE(0.5)*D_SQRT3*alpha[3]*beta[4]-DOUBLE(0.5)*D_SQRT3*alpha[6]*beta[5]+DOUBLE(0.5)*D_SQRT3*alpha[5]*beta[6];
            
            
        }
        // V.E.VDagger //
        void AdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *A,SU_Nc_ALGEBRA_FORMAT *ANew){
            
            ANew[0] = (3.0*A[2]*(imag(V[0])*imag(V[1]) - imag(V[3])*imag(V[4]) + real(V[0])*real(V[1]) - real(V[3])*real(V[4])) + D_SQRT3*A[7]*(imag(V[0])*imag(V[1]) + imag(V[3])*imag(V[4]) - 2.0*imag(V[6])*imag(V[7]) + real(V[0])*real(V[1]) + real(V[3])*real(V[4]) - 2.0*real(V[6])*real(V[7])) + 3.0*(A[5]*imag(V[4])*imag(V[6]) + A[5]*imag(V[3])*imag(V[7]) - A[1]*imag(V[4])*real(V[0]) - A[4]*imag(V[7])*real(V[0]) - A[1]*imag(V[3])*real(V[1]) - A[4]*imag(V[6])*real(V[1]) + A[1]*imag(V[1])*real(V[3]) - A[6]*imag(V[7])*real(V[3]) + A[1]*imag(V[0])*real(V[4]) - A[6]*imag(V[6])*real(V[4]) + A[0]*(imag(V[1])*imag(V[3]) + imag(V[0])*imag(V[4]) + real(V[1])*real(V[3]) + real(V[0])*real(V[4])) + A[4]*imag(V[1])*real(V[6]) + A[6]*imag(V[4])*real(V[6]) + A[5]*real(V[4])*real(V[6]) + A[4]*imag(V[0])*real(V[7]) + A[6]*imag(V[3])*real(V[7]) + A[5]*real(V[3])*real(V[7]) + A[3]*(imag(V[1])*imag(V[6]) + imag(V[0])*imag(V[7]) + real(V[1])*real(V[6]) + real(V[0])*real(V[7]))))/3.0;
            
            ANew[1] = -(A[6]*imag(V[4])*imag(V[6])) + A[6]*imag(V[3])*imag(V[7]) + A[2]*imag(V[1])*real(V[0]) + (A[7]*imag(V[1])*real(V[0]))/D_SQRT3 + A[0]*imag(V[4])*real(V[0]) + A[3]*imag(V[7])*real(V[0]) - A[2]*imag(V[0])*real(V[1]) - (A[7]*imag(V[0])*real(V[1]))/D_SQRT3 - A[0]*imag(V[3])*real(V[1]) - A[3]*imag(V[6])*real(V[1]) + A[0]*imag(V[1])*real(V[3]) - A[2]*imag(V[4])*real(V[3]) + (A[7]*imag(V[4])*real(V[3]))/D_SQRT3 + A[5]*imag(V[7])*real(V[3]) - A[0]*imag(V[0])*real(V[4]) + A[2]*imag(V[3])*real(V[4]) - (A[7]*imag(V[3])*real(V[4]))/D_SQRT3 - A[5]*imag(V[6])*real(V[4]) + A[1]*(-(imag(V[1])*imag(V[3])) + imag(V[0])*imag(V[4]) - real(V[1])*real(V[3]) + real(V[0])*real(V[4])) + A[3]*imag(V[1])*real(V[6]) + A[5]*imag(V[4])*real(V[6]) - (2.0*A[7]*imag(V[7])*real(V[6]))/D_SQRT3 - A[6]*real(V[4])*real(V[6]) - A[3]*imag(V[0])*real(V[7]) - A[5]*imag(V[3])*real(V[7]) + (2.0*A[7]*imag(V[6])*real(V[7]))/D_SQRT3 + A[6]*real(V[3])*real(V[7]) + A[4]*(-(imag(V[1])*imag(V[6])) + imag(V[0])*imag(V[7]) - real(V[1])*real(V[6]) + real(V[0])*real(V[7]));
            
            ANew[2] = (3.0*A[2]*(std::pow(imag(V[0]),2) - std::pow(imag(V[1]),2) - std::pow(imag(V[3]),2) + std::pow(imag(V[4]),2) + std::pow(real(V[0]),2) - std::pow(real(V[1]),2) - std::pow(real(V[3]),2) + std::pow(real(V[4]),2)) + D_SQRT3*A[7]* (std::pow(imag(V[0]),2) - std::pow(imag(V[1]),2) + std::pow(imag(V[3]),2) - std::pow(imag(V[4]),2) - 2.0*std::pow(imag(V[6]),2) + 2.0*std::pow(imag(V[7]),2) + std::pow(real(V[0]),2) - std::pow(real(V[1]),2) + std::pow(real(V[3]),2) - std::pow(real(V[4]),2) - 2.0*std::pow(real(V[6]),2) + 2.0*std::pow(real(V[7]),2)) + 6.0*(A[5]*imag(V[3])*imag(V[6]) - A[5]*imag(V[4])*imag(V[7]) - A[1]*imag(V[3])*real(V[0]) - A[4]*imag(V[6])*real(V[0]) + A[1]*imag(V[4])*real(V[1]) + A[4]*imag(V[7])*real(V[1]) + A[1]*imag(V[0])*real(V[3]) - A[6]*imag(V[6])*real(V[3]) - A[1]*imag(V[1])*real(V[4]) + A[6]*imag(V[7])*real(V[4]) + A[0]*(imag(V[0])*imag(V[3]) - imag(V[1])*imag(V[4]) + real(V[0])*real(V[3]) - real(V[1])*real(V[4])) + A[4]*imag(V[0])*real(V[6]) + A[6]*imag(V[3])*real(V[6]) + A[5]*real(V[3])*real(V[6]) - A[4]*imag(V[1])*real(V[7]) - A[6]*imag(V[4])*real(V[7]) - A[5]*real(V[4])*real(V[7]) + A[3]*(imag(V[0])*imag(V[6]) - imag(V[1])*imag(V[7]) + real(V[0])*real(V[6]) - real(V[1])*real(V[7]))))/6.0;
            
            ANew[3] = (3.0*A[2]*(imag(V[0])*imag(V[2]) - imag(V[3])*imag(V[5]) + real(V[0])*real(V[2]) - real(V[3])*real(V[5])) + D_SQRT3*A[7]*(imag(V[0])*imag(V[2]) + imag(V[3])*imag(V[5]) - 2.0*imag(V[6])*imag(V[8]) + real(V[0])*real(V[2]) + real(V[3])*real(V[5]) - 2.0*real(V[6])*real(V[8])) + 3.0*(A[5]*imag(V[5])*imag(V[6]) + A[5]*imag(V[3])*imag(V[8]) - A[1]*imag(V[5])*real(V[0]) - A[4]*imag(V[8])*real(V[0]) - A[1]*imag(V[3])*real(V[2]) - A[4]*imag(V[6])*real(V[2]) + A[1]*imag(V[2])*real(V[3]) - A[6]*imag(V[8])*real(V[3]) + A[1]*imag(V[0])*real(V[5]) - A[6]*imag(V[6])*real(V[5]) + A[0]*(imag(V[2])*imag(V[3]) + imag(V[0])*imag(V[5]) + real(V[2])*real(V[3]) + real(V[0])*real(V[5])) + A[4]*imag(V[2])*real(V[6]) + A[6]*imag(V[5])*real(V[6]) + A[5]*real(V[5])*real(V[6]) + A[4]*imag(V[0])*real(V[8]) + A[6]*imag(V[3])*real(V[8]) + A[5]*real(V[3])*real(V[8]) + A[3]*(imag(V[2])*imag(V[6]) + imag(V[0])*imag(V[8]) + real(V[2])*real(V[6]) + real(V[0])*real(V[8]))))/3.0;
            
            ANew[4] = -(A[6]*imag(V[5])*imag(V[6])) + A[6]*imag(V[3])*imag(V[8]) + A[2]*imag(V[2])*real(V[0]) + (A[7]*imag(V[2])*real(V[0]))/D_SQRT3 + A[0]*imag(V[5])*real(V[0]) + A[3]*imag(V[8])*real(V[0]) - A[2]*imag(V[0])*real(V[2]) - (A[7]*imag(V[0])*real(V[2]))/D_SQRT3 - A[0]*imag(V[3])*real(V[2]) - A[3]*imag(V[6])*real(V[2]) + A[0]*imag(V[2])*real(V[3]) - A[2]*imag(V[5])*real(V[3]) + (A[7]*imag(V[5])*real(V[3]))/D_SQRT3 + A[5]*imag(V[8])*real(V[3]) - A[0]*imag(V[0])*real(V[5]) + A[2]*imag(V[3])*real(V[5]) - (A[7]*imag(V[3])*real(V[5]))/D_SQRT3 - A[5]*imag(V[6])*real(V[5]) + A[1]*(-(imag(V[2])*imag(V[3])) + imag(V[0])*imag(V[5]) - real(V[2])*real(V[3]) + real(V[0])*real(V[5])) + A[3]*imag(V[2])*real(V[6]) + A[5]*imag(V[5])*real(V[6]) - (2.0*A[7]*imag(V[8])*real(V[6]))/D_SQRT3 - A[6]*real(V[5])*real(V[6]) - A[3]*imag(V[0])*real(V[8]) - A[5]*imag(V[3])*real(V[8]) + (2.0*A[7]*imag(V[6])*real(V[8]))/D_SQRT3 + A[6]*real(V[3])*real(V[8]) + A[4]*(-(imag(V[2])*imag(V[6])) + imag(V[0])*imag(V[8]) - real(V[2])*real(V[6]) + real(V[0])*real(V[8]));
            
            ANew[5] = (3.0*A[2]*(imag(V[1])*imag(V[2]) - imag(V[4])*imag(V[5]) + real(V[1])*real(V[2]) - real(V[4])*real(V[5])) + D_SQRT3*A[7]*(imag(V[1])*imag(V[2]) + imag(V[4])*imag(V[5]) - 2.0*imag(V[7])*imag(V[8]) + real(V[1])*real(V[2]) + real(V[4])*real(V[5]) - 2.0*real(V[7])*real(V[8])) + 3.0*(A[5]*imag(V[5])*imag(V[7]) + A[5]*imag(V[4])*imag(V[8]) - A[1]*imag(V[5])*real(V[1]) - A[4]*imag(V[8])*real(V[1]) - A[1]*imag(V[4])*real(V[2]) - A[4]*imag(V[7])*real(V[2]) + A[1]*imag(V[2])*real(V[4]) - A[6]*imag(V[8])*real(V[4]) + A[1]*imag(V[1])*real(V[5]) - A[6]*imag(V[7])*real(V[5]) + A[0]*(imag(V[2])*imag(V[4]) + imag(V[1])*imag(V[5]) + real(V[2])*real(V[4]) + real(V[1])*real(V[5])) + A[4]*imag(V[2])*real(V[7]) + A[6]*imag(V[5])*real(V[7]) + A[5]*real(V[5])*real(V[7]) + A[4]*imag(V[1])*real(V[8]) + A[6]*imag(V[4])*real(V[8]) + A[5]*real(V[4])*real(V[8]) + A[3]*(imag(V[2])*imag(V[7]) + imag(V[1])*imag(V[8]) + real(V[2])*real(V[7]) + real(V[1])*real(V[8]))))/3.0;
            
            ANew[6] = -(A[6]*imag(V[5])*imag(V[7])) + A[6]*imag(V[4])*imag(V[8]) + A[2]*imag(V[2])*real(V[1]) + (A[7]*imag(V[2])*real(V[1]))/D_SQRT3 + A[0]*imag(V[5])*real(V[1]) + A[3]*imag(V[8])*real(V[1]) - A[2]*imag(V[1])*real(V[2]) - (A[7]*imag(V[1])*real(V[2]))/D_SQRT3 - A[0]*imag(V[4])*real(V[2]) - A[3]*imag(V[7])*real(V[2]) + A[0]*imag(V[2])*real(V[4]) - A[2]*imag(V[5])*real(V[4]) + (A[7]*imag(V[5])*real(V[4]))/D_SQRT3 + A[5]*imag(V[8])*real(V[4]) - A[0]*imag(V[1])*real(V[5]) + A[2]*imag(V[4])*real(V[5]) - (A[7]*imag(V[4])*real(V[5]))/D_SQRT3 - A[5]*imag(V[7])*real(V[5]) + A[1]*(-(imag(V[2])*imag(V[4])) + imag(V[1])*imag(V[5]) - real(V[2])*real(V[4]) + real(V[1])*real(V[5])) + A[3]*imag(V[2])*real(V[7]) + A[5]*imag(V[5])*real(V[7]) - (2.0*A[7]*imag(V[8])*real(V[7]))/D_SQRT3 - A[6]*real(V[5])*real(V[7]) - A[3]*imag(V[1])*real(V[8]) - A[5]*imag(V[4])*real(V[8]) + (2.0*A[7]*imag(V[7])*real(V[8]))/D_SQRT3 + A[6]*real(V[4])*real(V[8]) + A[4]*(-(imag(V[2])*imag(V[7])) + imag(V[1])*imag(V[8]) - real(V[2])*real(V[7]) + real(V[1])*real(V[8]));
            
            ANew[7] = (D_SQRT3*A[2]*(std::pow(imag(V[0]),2) + std::pow(imag(V[1]),2) - 2.0*std::pow(imag(V[2]),2) - std::pow(imag(V[3]),2) - std::pow(imag(V[4]),2) + 2.0*std::pow(imag(V[5]),2) +std::pow(real(V[0]),2) + std::pow(real(V[1]),2) - 2.0*std::pow(real(V[2]),2) - std::pow(real(V[3]),2) - std::pow(real(V[4]),2) + 2.0*std::pow(real(V[5]),2)) + A[7]*(std::pow(imag(V[0]),2) + std::pow(imag(V[1]),2) - 2.0*std::pow(imag(V[2]),2) + std::pow(imag(V[3]),2) + std::pow(imag(V[4]),2) - 2.0*std::pow(imag(V[5]),2) - 2.0*std::pow(imag(V[6]),2) -2.0*std::pow(imag(V[7]),2) + 4*std::pow(imag(V[8]),2) + std::pow(real(V[0]),2) + std::pow(real(V[1]),2) - 2.0*std::pow(real(V[2]),2) + std::pow(real(V[3]),2) + std::pow(real(V[4]),2) - 2.0*std::pow(real(V[5]),2) - 2.0*std::pow(real(V[6]),2) - 2.0*std::pow(real(V[7]),2) + 4*std::pow(real(V[8]),2)) + 2.0*D_SQRT3*(A[5]*imag(V[3])*imag(V[6]) + A[5]*imag(V[4])*imag(V[7]) - 2.0*A[5]*imag(V[5])*imag(V[8]) - A[1]*imag(V[3])*real(V[0]) - A[4]*imag(V[6])*real(V[0]) - A[1]*imag(V[4])*real(V[1]) - A[4]*imag(V[7])*real(V[1]) + 2.0*A[1]*imag(V[5])*real(V[2]) + 2.0*A[4]*imag(V[8])*real(V[2]) + A[1]*imag(V[0])*real(V[3]) - A[6]*imag(V[6])*real(V[3]) + A[1]*imag(V[1])*real(V[4]) - A[6]*imag(V[7])*real(V[4]) - 2.0*A[1]*imag(V[2])*real(V[5]) + 2.0*A[6]*imag(V[8])*real(V[5]) + A[0]*(imag(V[0])*imag(V[3]) + imag(V[1])*imag(V[4]) - 2.0*imag(V[2])*imag(V[5]) + real(V[0])*real(V[3]) + real(V[1])*real(V[4]) - 2.0*real(V[2])*real(V[5])) + A[4]*imag(V[0])*real(V[6]) + A[6]*imag(V[3])*real(V[6]) + A[5]*real(V[3])*real(V[6]) + A[4]*imag(V[1])*real(V[7]) + A[6]*imag(V[4])*real(V[7]) + A[5]*real(V[4])*real(V[7]) - 2.0*A[4]*imag(V[2])*real(V[8]) - 2.0*A[6]*imag(V[5])*real(V[8]) - 2.0*A[5]*real(V[5])*real(V[8]) + A[3]*(imag(V[0])*imag(V[6]) + imag(V[1])*imag(V[7]) - 2.0*imag(V[2])*imag(V[8]) + real(V[0])*real(V[6]) + real(V[1])*real(V[7]) - 2.0*real(V[2])*real(V[8]))))/6.0;
            
        }
        
        // VDagger.E.V //
        void InverseAdjointMultiplication(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ALGEBRA_FORMAT *A,SU_Nc_ALGEBRA_FORMAT *ANew){
            
            ANew[0] = (3.0*A[2]*(imag(V[0])*imag(V[3]) - imag(V[1])*imag(V[4]) + real(V[0])*real(V[3]) - real(V[1])*real(V[4])) + D_SQRT3*A[7]*(imag(V[0])*imag(V[3]) + imag(V[1])*imag(V[4]) - 2.0*imag(V[2])*imag(V[5]) + real(V[0])*real(V[3]) + real(V[1])*real(V[4]) - 2.0*real(V[2])*real(V[5])) + 3.0*(A[5]*imag(V[2])*imag(V[4]) + A[5]*imag(V[1])*imag(V[5]) + A[1]*imag(V[4])*real(V[0]) + A[4]*imag(V[5])*real(V[0]) - A[1]*imag(V[3])*real(V[1]) + A[6]*imag(V[5])*real(V[1]) - A[4]*imag(V[3])*real(V[2]) - A[6]*imag(V[4])*real(V[2]) + A[1]*imag(V[1])*real(V[3]) + A[4]*imag(V[2])*real(V[3]) - A[1]*imag(V[0])*real(V[4]) + A[6]*imag(V[2])*real(V[4]) + A[5]*real(V[2])*real(V[4]) + A[0]*(imag(V[1])*imag(V[3]) + imag(V[0])*imag(V[4]) + real(V[1])*real(V[3]) + real(V[0])*real(V[4])) - A[4]*imag(V[0])*real(V[5]) - A[6]*imag(V[1])*real(V[5]) + A[5]*real(V[1])*real(V[5]) + A[3]*(imag(V[2])*imag(V[3]) + imag(V[0])*imag(V[5]) + real(V[2])*real(V[3]) + real(V[0])*real(V[5]))))/3.0;
            
            ANew[1] = -(A[6]*imag(V[2])*imag(V[4])) + A[6]*imag(V[1])*imag(V[5]) - A[2]*imag(V[3])*real(V[0]) - (A[7]*imag(V[3])*real(V[0]))/D_SQRT3 - A[0]*imag(V[4])*real(V[0]) - A[3]*imag(V[5])*real(V[0]) - A[0]*imag(V[3])*real(V[1]) + A[2]*imag(V[4])*real(V[1]) - (A[7]*imag(V[4])*real(V[1]))/D_SQRT3 - A[5]*imag(V[5])*real(V[1]) - A[3]*imag(V[3])*real(V[2]) - A[5]*imag(V[4])*real(V[2]) + (2.0*A[7]*imag(V[5])*real(V[2]))/D_SQRT3 + A[2]*imag(V[0])*real(V[3]) + (A[7]*imag(V[0])*real(V[3]))/D_SQRT3 + A[0]*imag(V[1])*real(V[3]) + A[3]*imag(V[2])*real(V[3]) + A[0]*imag(V[0])*real(V[4]) - A[2]*imag(V[1])*real(V[4]) + (A[7]*imag(V[1])*real(V[4]))/D_SQRT3 + A[5]*imag(V[2])*real(V[4]) - A[6]*real(V[2])*real(V[4]) + A[1]*(-(imag(V[1])*imag(V[3])) + imag(V[0])*imag(V[4]) - real(V[1])*real(V[3]) + real(V[0])*real(V[4])) + A[3]*imag(V[0])*real(V[5]) + A[5]*imag(V[1])*real(V[5]) - (2.0*A[7]*imag(V[2])*real(V[5]))/D_SQRT3 + A[6]*real(V[1])*real(V[5]) + A[4]*(-(imag(V[2])*imag(V[3])) + imag(V[0])*imag(V[5]) - real(V[2])*real(V[3]) + real(V[0])*real(V[5]));
            
            ANew[2] = (3.0*A[2]*(std::pow(imag(V[0]),2) - std::pow(imag(V[1]),2) - std::pow(imag(V[3]),2) + std::pow(imag(V[4]),2) + std::pow(real(V[0]),2) - std::pow(real(V[1]),2) - std::pow(real(V[3]),2) + std::pow(real(V[4]),2)) + D_SQRT3*A[7]* (std::pow(imag(V[0]),2) + std::pow(imag(V[1]),2) - 2.0*std::pow(imag(V[2]),2) - std::pow(imag(V[3]),2) - std::pow(imag(V[4]),2) + 2.0*std::pow(imag(V[5]),2) + std::pow(real(V[0]),2) + std::pow(real(V[1]),2) - 2.0*std::pow(real(V[2]),2) - std::pow(real(V[3]),2) - std::pow(real(V[4]),2) + 2.0*std::pow(real(V[5]),2)) + 6.0*(A[5]*imag(V[1])*imag(V[2]) - A[5]*imag(V[4])*imag(V[5]) + A[1]*imag(V[1])*real(V[0]) + A[4]*imag(V[2])*real(V[0]) - A[1]*imag(V[0])*real(V[1]) + A[6]*imag(V[2])*real(V[1]) - A[4]*imag(V[0])*real(V[2]) - A[6]*imag(V[1])*real(V[2]) + A[5]*real(V[1])*real(V[2]) - A[1]*imag(V[4])*real(V[3]) - A[4]*imag(V[5])*real(V[3]) + A[1]*imag(V[3])*real(V[4]) - A[6]*imag(V[5])*real(V[4]) + A[0]*(imag(V[0])*imag(V[1]) - imag(V[3])*imag(V[4]) + real(V[0])*real(V[1]) - real(V[3])*real(V[4])) + A[4]*imag(V[3])*real(V[5]) + A[6]*imag(V[4])*real(V[5]) - A[5]*real(V[4])*real(V[5]) + A[3]*(imag(V[0])*imag(V[2]) - imag(V[3])*imag(V[5]) + real(V[0])*real(V[2]) - real(V[3])*real(V[5]))))/6.0;
            
            ANew[3] = (3.0*A[2]*(imag(V[0])*imag(V[6]) - imag(V[1])*imag(V[7]) + real(V[0])*real(V[6]) - real(V[1])*real(V[7])) + D_SQRT3*A[7]*(imag(V[0])*imag(V[6]) + imag(V[1])*imag(V[7]) - 2.0*imag(V[2])*imag(V[8]) + real(V[0])*real(V[6]) + real(V[1])*real(V[7]) - 2.0*real(V[2])*real(V[8])) + 3.0*(A[5]*imag(V[2])*imag(V[7]) + A[5]*imag(V[1])*imag(V[8]) + A[1]*imag(V[7])*real(V[0]) + A[4]*imag(V[8])*real(V[0]) - A[1]*imag(V[6])*real(V[1]) + A[6]*imag(V[8])*real(V[1]) - A[4]*imag(V[6])*real(V[2]) - A[6]*imag(V[7])*real(V[2]) + A[1]*imag(V[1])*real(V[6]) + A[4]*imag(V[2])*real(V[6]) - A[1]*imag(V[0])*real(V[7]) + A[6]*imag(V[2])*real(V[7]) + A[5]*real(V[2])*real(V[7]) + A[0]*(imag(V[1])*imag(V[6]) + imag(V[0])*imag(V[7]) + real(V[1])*real(V[6]) + real(V[0])*real(V[7])) - A[4]*imag(V[0])*real(V[8]) - A[6]*imag(V[1])*real(V[8]) + A[5]*real(V[1])*real(V[8]) + A[3]*(imag(V[2])*imag(V[6]) + imag(V[0])*imag(V[8]) + real(V[2])*real(V[6]) + real(V[0])*real(V[8]))))/3.0;
            
            ANew[4] = -(A[6]*imag(V[2])*imag(V[7])) + A[6]*imag(V[1])*imag(V[8]) - A[2]*imag(V[6])*real(V[0]) - (A[7]*imag(V[6])*real(V[0]))/D_SQRT3 - A[0]*imag(V[7])*real(V[0]) - A[3]*imag(V[8])*real(V[0]) - A[0]*imag(V[6])*real(V[1]) + A[2]*imag(V[7])*real(V[1]) - (A[7]*imag(V[7])*real(V[1]))/D_SQRT3 - A[5]*imag(V[8])*real(V[1]) - A[3]*imag(V[6])*real(V[2]) - A[5]*imag(V[7])*real(V[2]) + (2.0*A[7]*imag(V[8])*real(V[2]))/D_SQRT3 + A[2]*imag(V[0])*real(V[6]) + (A[7]*imag(V[0])*real(V[6]))/D_SQRT3 + A[0]*imag(V[1])*real(V[6]) + A[3]*imag(V[2])*real(V[6]) + A[0]*imag(V[0])*real(V[7]) - A[2]*imag(V[1])*real(V[7]) + (A[7]*imag(V[1])*real(V[7]))/D_SQRT3 + A[5]*imag(V[2])*real(V[7]) - A[6]*real(V[2])*real(V[7]) + A[1]*(-(imag(V[1])*imag(V[6])) + imag(V[0])*imag(V[7]) - real(V[1])*real(V[6]) + real(V[0])*real(V[7])) + A[3]*imag(V[0])*real(V[8]) + A[5]*imag(V[1])*real(V[8]) - (2.0*A[7]*imag(V[2])*real(V[8]))/D_SQRT3 + A[6]*real(V[1])*real(V[8]) + A[4]*(-(imag(V[2])*imag(V[6])) + imag(V[0])*imag(V[8]) - real(V[2])*real(V[6]) + real(V[0])*real(V[8]));
            
            ANew[5] = (3.0*A[2]*(imag(V[3])*imag(V[6]) - imag(V[4])*imag(V[7]) + real(V[3])*real(V[6]) - real(V[4])*real(V[7])) + D_SQRT3*A[7]*(imag(V[3])*imag(V[6]) + imag(V[4])*imag(V[7]) - 2.0*imag(V[5])*imag(V[8]) + real(V[3])*real(V[6]) + real(V[4])*real(V[7]) - 2.0*real(V[5])*real(V[8])) + 3.0*(A[5]*imag(V[5])*imag(V[7]) + A[5]*imag(V[4])*imag(V[8]) + A[1]*imag(V[7])*real(V[3]) + A[4]*imag(V[8])*real(V[3]) - A[1]*imag(V[6])*real(V[4]) + A[6]*imag(V[8])*real(V[4]) - A[4]*imag(V[6])*real(V[5]) - A[6]*imag(V[7])*real(V[5]) + A[1]*imag(V[4])*real(V[6]) + A[4]*imag(V[5])*real(V[6]) - A[1]*imag(V[3])*real(V[7]) + A[6]*imag(V[5])*real(V[7]) + A[5]*real(V[5])*real(V[7]) + A[0]*(imag(V[4])*imag(V[6]) + imag(V[3])*imag(V[7]) + real(V[4])*real(V[6]) + real(V[3])*real(V[7])) - A[4]*imag(V[3])*real(V[8]) - A[6]*imag(V[4])*real(V[8]) + A[5]*real(V[4])*real(V[8]) + A[3]*(imag(V[5])*imag(V[6]) + imag(V[3])*imag(V[8]) + real(V[5])*real(V[6]) + real(V[3])*real(V[8]))))/3.0;
            
            ANew[6] = -(A[6]*imag(V[5])*imag(V[7])) + A[6]*imag(V[4])*imag(V[8]) - A[2]*imag(V[6])*real(V[3]) - (A[7]*imag(V[6])*real(V[3]))/D_SQRT3 - A[0]*imag(V[7])*real(V[3]) - A[3]*imag(V[8])*real(V[3]) - A[0]*imag(V[6])*real(V[4]) + A[2]*imag(V[7])*real(V[4]) - (A[7]*imag(V[7])*real(V[4]))/D_SQRT3 - A[5]*imag(V[8])*real(V[4]) - A[3]*imag(V[6])*real(V[5]) - A[5]*imag(V[7])*real(V[5]) + (2.0*A[7]*imag(V[8])*real(V[5]))/D_SQRT3 + A[2]*imag(V[3])*real(V[6]) + (A[7]*imag(V[3])*real(V[6]))/D_SQRT3 + A[0]*imag(V[4])*real(V[6]) + A[3]*imag(V[5])*real(V[6]) + A[0]*imag(V[3])*real(V[7]) - A[2]*imag(V[4])*real(V[7]) + (A[7]*imag(V[4])*real(V[7]))/D_SQRT3 + A[5]*imag(V[5])*real(V[7]) - A[6]*real(V[5])*real(V[7]) + A[1]*(-(imag(V[4])*imag(V[6])) + imag(V[3])*imag(V[7]) - real(V[4])*real(V[6]) + real(V[3])*real(V[7])) + A[3]*imag(V[3])*real(V[8]) + A[5]*imag(V[4])*real(V[8]) - (2.0*A[7]*imag(V[5])*real(V[8]))/D_SQRT3 + A[6]*real(V[4])*real(V[8]) + A[4]*(-(imag(V[5])*imag(V[6])) + imag(V[3])*imag(V[8]) - real(V[5])*real(V[6]) + real(V[3])*real(V[8]));
            
            ANew[7] = (D_SQRT3*A[2]*(std::pow(imag(V[0]),2) - std::pow(imag(V[1]),2) + std::pow(imag(V[3]),2) - std::pow(imag(V[4]),2) - 2.0*std::pow(imag(V[6]),2) + 2.0*std::pow(imag(V[7]),2) + std::pow(real(V[0]),2) - std::pow(real(V[1]),2) + std::pow(real(V[3]),2) - std::pow(real(V[4]),2) - 2.0*std::pow(real(V[6]),2) + 2.0*std::pow(real(V[7]),2)) + A[7]*(std::pow(imag(V[0]),2) + std::pow(imag(V[1]),2) - 2.0*std::pow(imag(V[2]),2) + std::pow(imag(V[3]),2) + std::pow(imag(V[4]),2) - 2.0*std::pow(imag(V[5]),2) - 2.0*std::pow(imag(V[6]),2) - 2.0*std::pow(imag(V[7]),2) + 4.0*std::pow(imag(V[8]),2) + std::pow(real(V[0]),2) + std::pow(real(V[1]),2) - 2.0*std::pow(real(V[2]),2) + std::pow(real(V[3]),2) + std::pow(real(V[4]),2) - 2.0*std::pow(real(V[5]),2) - 2.0*std::pow(real(V[6]),2) - 2.0*std::pow(real(V[7]),2) + 4.0*std::pow(real(V[8]),2)) + 2.0*D_SQRT3*(A[5]*imag(V[1])*imag(V[2]) + A[5]*imag(V[4])*imag(V[5]) - 2.0*A[5]*imag(V[7])*imag(V[8]) + A[1]*imag(V[1])*real(V[0]) + A[4]*imag(V[2])*real(V[0]) - A[1]*imag(V[0])*real(V[1]) + A[6]*imag(V[2])*real(V[1]) - A[4]*imag(V[0])*real(V[2]) - A[6]*imag(V[1])*real(V[2]) + A[5]*real(V[1])*real(V[2]) + A[1]*imag(V[4])*real(V[3]) + A[4]*imag(V[5])*real(V[3]) - A[1]*imag(V[3])*real(V[4]) + A[6]*imag(V[5])*real(V[4]) - A[4]*imag(V[3])*real(V[5]) - A[6]*imag(V[4])*real(V[5]) + A[5]*real(V[4])*real(V[5]) - 2.0*A[1]*imag(V[7])*real(V[6]) - 2.0*A[4]*imag(V[8])*real(V[6]) + 2.0*A[1]*imag(V[6])*real(V[7]) - 2.0*A[6]*imag(V[8])*real(V[7]) + A[0]*(imag(V[0])*imag(V[1]) + imag(V[3])*imag(V[4]) - 2.0*imag(V[6])*imag(V[7]) + real(V[0])*real(V[1]) + real(V[3])*real(V[4]) - 2.0*real(V[6])*real(V[7])) + 2.0*A[4]*imag(V[6])*real(V[8]) + 2.0*A[6]*imag(V[7])*real(V[8]) - 2.0*A[5]*real(V[7])*real(V[8]) + A[3]* (imag(V[0])*imag(V[2]) + imag(V[3])*imag(V[5]) - 2.0*imag(V[6])*imag(V[8]) + real(V[0])*real(V[2]) + real(V[3])*real(V[5]) - 2.0*real(V[6])*real(V[8]))))/6.0;
        }
        
        //MATRIX TO STRING REPRESENTATION //
        std::string VectorToString(SU_Nc_ALGEBRA_FORMAT *A){
            
            std::stringstream sstm;
            
            sstm.precision(OUTPUT_PRECISION);
            
            for(int alpha=0;alpha<VectorSize;alpha++){
                sstm << A[alpha] << " ";
            }
            
            return sstm.str();
            
        }

    }
    
}

#endif
