#ifndef __SU_Nc__GROUP_OPERATIONS__
#define __SU_Nc__GROUP_OPERATIONS__

namespace SUNcGroup{
    
    //MATRIX SIZE
    static const int MatrixSize=1;
    
    //UNIT MATRIX
    static const SU_Nc_FUNDAMENTAL_FORMAT UnitMatrix[MatrixSize]={DOUBLE(1.0)};
    
    //INDEX CONVENTIONS
    int AdjIndex(int a,int b){
        return a+1*b;
    }
    
    //OPERATIONS INVOLVING U(1) MATRICES
    namespace Operations{
        
        COMPLEX Det(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return V[0];
        }
        
        /////////////////////
        //INVERSE          //
        /////////////////////
        
        void Inverse(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VDagger){
            VDagger[0]=conj(V[0]);
        }
        
        /////////////////////
        //TRACE OPERATIONS //
        /////////////////////
        
        //tr(V)
        COMPLEX tr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return V[0];
        }
        
        DOUBLE ReTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return real(V[0]);
        }
        
        DOUBLE ImTr(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return imag(V[0]);
        }
        
        
        //tr(1-V)
        COMPLEX trIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(1.0)-V[0];
        }
        
        DOUBLE ReTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(1.0)-real(V[0]);
        }
        
        DOUBLE ImTrIDMinusU(SU_Nc_FUNDAMENTAL_FORMAT *V){
            return DOUBLE(1.0)-imag(V[0]);
        }
        
        //tr(it^{a}V)
        void trIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,COMPLEX *trItV){
            
            trItV[0]= DOUBLE(0.5)*M_SQRT2*ComplexI*V[0];
            
        }
        
        void ReTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-DOUBLE(0.5)*M_SQRT2*imag(V[0]);
            
        }
        
        void ReTrIGenU(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *RetrItV){
            
            RetrItV[0]=-c*DOUBLE(0.5)*M_SQRT2*imag(V[0]);
            
        }
        
        void ImTrIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE *ImtrItV){
            
            ImtrItV[0]=DOUBLE(0.5)*M_SQRT2*real(V[0]);
            
        }
        
        
        
        //BASIC MATRIX MULTIPLICATIONS
        void UU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2){
            
            V1V2[0]=V1[0]*V2[0];
            
        }
        
        
        void UD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1V2D){
            
            V1V2D[0]=V1[0]*conj(V2[0]);
            
        }
        
        void DU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2){
            
            V1DV2[0]=conj(V1[0])*V2[0];
            
        }
        
        void DD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2D){
            
            V1DV2D[0]=conj(V1[0]*V2[0]);
            
        }
        
        
        //CONSTRUCTION OF ADJOINT REPRESENTATION U^{adj}_{ab}=2 tr[t^{a} U_{f} t^{b} U_{f}^{\dagger}]  
        void GetAdjoint(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_ADJOINT_FORMAT *VAdj){
            VAdj[0]=1.0;
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
        
        void ReTrIGenIGenU(SU_Nc_FUNDAMENTAL_FORMAT *V,DOUBLE ReTrItItV[1][1]){
            ReTrItItV[0][0]=-DOUBLE(0.5)*real(V[0]);
        }
        
    }
    
    namespace Extended{
        
        //COMPUTE UNITARIZATION UTilde=U/Sqrt(U^{dagger} U) //
        void MaxTraceProjection(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VTilde){
            VTilde[0]=(V[0]/abs(V[0]));
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
