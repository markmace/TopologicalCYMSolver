#ifndef __SU_Nc__GROUP_ADVANCED_OPERATIONS__
#define __SU_Nc__GROUP_ADVANCED_OPERATIONS__

#include <string>
#include <sstream>

namespace SUNcGroup{ 
    
    ////////////////////////////////////
    //MATRIX SUM                      //
    ////////////////////////////////////
    
    
    void MatrixSum(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1PlusV2){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            V1PlusV2[s]=V1[s]+V2[s];
        }
    }
    
    void MatrixSum(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,const SU_Nc_FUNDAMENTAL_FORMAT *V3,const SU_Nc_FUNDAMENTAL_FORMAT *V4, SU_Nc_FUNDAMENTAL_FORMAT *Sum){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            Sum[s]=(V1[s]+V2[s]+V3[s]+V4[s]);
        }
    }
    
    void MatrixSum(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,const SU_Nc_FUNDAMENTAL_FORMAT *V3,const SU_Nc_FUNDAMENTAL_FORMAT *V4,const SU_Nc_FUNDAMENTAL_FORMAT *V5,const SU_Nc_FUNDAMENTAL_FORMAT *V6,const SU_Nc_FUNDAMENTAL_FORMAT *V7,const SU_Nc_FUNDAMENTAL_FORMAT *V8, SU_Nc_FUNDAMENTAL_FORMAT *Sum){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            Sum[s]=(V1[s]+V2[s]+V3[s]+V4[s]+V5[s]+V6[s]+V7[s]+V8[s]);
        }
    }
    
    void AvgMatrix(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,const SU_Nc_FUNDAMENTAL_FORMAT *V3,const SU_Nc_FUNDAMENTAL_FORMAT *V4, SU_Nc_FUNDAMENTAL_FORMAT *Sum){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            Sum[s]=DOUBLE(0.25)*(V1[s]+V2[s]+V3[s]+V4[s]);
        }
    }
    
    void AvgMatrix(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,const SU_Nc_FUNDAMENTAL_FORMAT *V3,const SU_Nc_FUNDAMENTAL_FORMAT *V4,const SU_Nc_FUNDAMENTAL_FORMAT *V5,const SU_Nc_FUNDAMENTAL_FORMAT *V6,const SU_Nc_FUNDAMENTAL_FORMAT *V7,const SU_Nc_FUNDAMENTAL_FORMAT *V8, SU_Nc_FUNDAMENTAL_FORMAT *Sum){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            Sum[s]=DOUBLE(0.125)*(V1[s]+V2[s]+V3[s]+V4[s]+V5[s]+V6[s]+V7[s]+V8[s]);
        }
    }
    
    void MatrixDifference(const SU_Nc_FUNDAMENTAL_FORMAT *V1,const SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V1MinusV2){
        
        for(int s=0;s<SUNcGroup::MatrixSize;s++){
            V1MinusV2[s]=V1[s]-V2[s];
        }
    }
    
    namespace AdvancedOperations{
        
        ////////////////////
        // POWER FUNCTION //
        ////////////////////
        
        void Power(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *Vm,DOUBLE m){
            
            SU_Nc_ALGEBRA_FORMAT logV[SUNcAlgebra::VectorSize];
            
            SUNcAlgebra::Operations::MatrixILog(DOUBLE(1.0),V,logV);
            
            SUNcAlgebra::Operations::MatrixIExp(m,logV,Vm);
            
        }
        
        
        ////////////////////////////////////
        //COMPOSITE MATRIX MULTIPLICATIONS//
        ////////////////////////////////////
        
        //V1.V2.V3D
        void UUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2,V3,V1V2V3);
            
        }
        
        //V1.V2.V3D 
        void UUD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; 
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UD(V1V2,V3,V1V2V3D);
            
        }
        
        //V1.V2D.V3D
        void UDD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V1V2DV3D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UD(V1,V2,V1V2D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UD(V1V2D,V3,V1V2DV3D);
            
        }
        
        //V1D.V2.V3 
        void DUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; 
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2,V3,V1DV2V3);
            
        }
        
        //V1D.V2D.V3 
        void DDU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2DV3){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2D[SUNcGroup::MatrixSize]; 
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DD(V1,V2,V1DV2D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2D,V3,V1DV2DV3);
            
        }
        
        //V1.V2.V3.V4
        void UUUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3V4){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::UU(V3,V4,V3V4);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2,V3V4,V1V2V3V4);
            
        }
        
        //V1.V2.V3.V4D
        void UUUD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3V4D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::UD(V3,V4,V3V4D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2,V3V4D,V1V2V3V4D);
            
        }
        
        //V1D.V2.V3.V4
        void DUUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3V4){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2); SUNcGroup::Operations::UU(V3,V4,V3V4);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2,V3V4,V1DV2V3V4);
            
        }
        
        //V1.V2.V3D.V4D    
        void UUDD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3DV4D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3DV4D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::DD(V3,V4,V3DV4D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2,V3DV4D,V1V2V3DV4D);
            
        }
        
        //V1D.V2D.V3.V4    
        void DDUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2DV3V4){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DD(V1,V2,V1DV2D); SUNcGroup::Operations::UU(V3,V4,V3V4);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2D,V3V4,V1DV2DV3V4);
            
        }
        
        //V1.V2D.V3D.V4    
        void UDDU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1V2DV3DV4){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3DV4[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UD(V1,V2,V1V2D); SUNcGroup::Operations::DU(V3,V4,V3DV4);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2D,V3DV4,V1V2DV3DV4);
            
        }
        
        //V1D.V2.V3.V4D    
        void DUUD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3V4D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2); SUNcGroup::Operations::UD(V3,V4,V3V4D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2,V3V4D,V1DV2V3V4D);
            
        }

        //V1.V2.V3.V4D.V5D.V6D
        void UUUDDD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3V4DV5DV6D){
            
            
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3V4D[SUNcGroup::MatrixSize];  SU_Nc_FUNDAMENTAL_FORMAT V5DV6D[SUNcGroup::MatrixSize];
            
            
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::UD(V3,V4,V3V4D); SUNcGroup::Operations::DD(V5,V6,V5DV6D);
            
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::AdvancedOperations::UUU(V1V2,V3V4D,V5DV6D,V1V2V3V4DV5DV6D);
        }
        
        //V1D.V2D.V3D.V4.V5.V6
        void DDDUUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2DV3DV4V5V6){
            
            
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V3DV4[SUNcGroup::MatrixSize];  SU_Nc_FUNDAMENTAL_FORMAT V5V6[SUNcGroup::MatrixSize];
            
            
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DD(V1,V2,V1DV2D); SUNcGroup::Operations::DU(V3,V4,V3DV4); SUNcGroup::Operations::UU(V5,V6,V5V6);
            
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::AdvancedOperations::UUU(V1DV2D,V3DV4,V5V6,V1DV2DV3DV4V5V6);
            
            
            
        }
        //V1.V2.V3D.V4D.V5D.V6
        void UUDDDU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3DV4DV5DV6){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5DV6[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1V2V3D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4DV5DV6[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::DU(V5,V6,V5DV6);
            SUNcGroup::Operations::UD(V1V2,V3,V1V2V3D); SUNcGroup::Operations::DU(V4,V5DV6,V4DV5DV6);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2V3D,V4DV5DV6,V1V2V3DV4DV5DV6);
        }
        //V1D.V2D.V3.V4.V5.V6D
        void DDUUUD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2DV3V4V5V6D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5V6D[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1DV2DV3[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4V5V6D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DD(V1,V2,V1DV2D); SUNcGroup::Operations::UD(V5,V6,V5V6D);
            SUNcGroup::Operations::UU(V1DV2D,V3,V1DV2DV3); SUNcGroup::Operations::UU(V4,V5V6D,V4V5V6D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2DV3,V4V5V6D,V1DV2DV3V4V5V6D);
            
        }
        //V1.V2.V3D.V4D.V5.V6D
        void UUDDUD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1V2V3DV4DV5V6D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5V6D[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1V2V3D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4DV5V6D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UU(V1,V2,V1V2); SUNcGroup::Operations::UD(V5,V6,V5V6D);
            SUNcGroup::Operations::UD(V1V2,V3,V1V2V3D); SUNcGroup::Operations::DU(V4,V5V6D,V4DV5V6D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2V3D,V4DV5V6D,V1V2V3DV4DV5V6D);
            
        }
        
        //V1.V2D.V3.V4.V5D.V6D
        void UDUUDD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1V2DV3V4V5DV6D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5DV6D[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1V2DV3[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4V5DV6D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UD(V1,V2,V1V2D); SUNcGroup::Operations::DD(V5,V6,V5DV6D);
            SUNcGroup::Operations::UU(V1V2D,V3,V1V2DV3); SUNcGroup::Operations::UU(V4,V5DV6D,V4V5DV6D);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2DV3,V4V5DV6D,V1V2DV3V4V5DV6D);
        }
        //V1D.V2.V3D.V4D.V5.V6
        void DUDDUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3DV4DV5V6){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5V6[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1DV2V3D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4DV5V6[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2); SUNcGroup::Operations::UU(V5,V6,V5V6);
            SUNcGroup::Operations::UD(V1DV2,V3,V1DV2V3D); SUNcGroup::Operations::DU(V4,V5V6,V4DV5V6);
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2V3D,V4DV5V6,V1DV2V3DV4DV5V6);
            
        }
        
        //V1.V2D.V3D.V4D.V5.V6
        void UDDDUU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1V2DV3DV4DV5V6){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1V2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5V6[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1V2DV3D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4DV5V6[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::UD(V1,V2,V1V2D); SUNcGroup::Operations::UU(V5,V6,V5V6);
            SUNcGroup::Operations::UD(V1V2D,V3,V1V2DV3D); SUNcGroup::Operations::DU(V4,V5V6,V4DV5V6);
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1V2DV3D,V4DV5V6,V1V2DV3DV4DV5V6);
        }
        //V1D.V2.V3.V4.V5D.V6D
        void DUUUDD(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3V4V5DV6D){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5DV6D[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1DV2V3[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4V5DV6D[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2); SUNcGroup::Operations::DD(V5,V6,V5DV6D);
            SUNcGroup::Operations::UU(V1DV2,V3,V1DV2V3); SUNcGroup::Operations::UU(V4,V5DV6D,V4V5DV6D);
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2V3,V4V5DV6D,V1DV2V3V4V5DV6D);
            
        }
        //V1D.V2D.V3.V4.V5D.V6
        void DDUUDU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2DV3V4V5DV6){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2D[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5DV6[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1DV2DV3[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4V5DV6[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DD(V1,V2,V1DV2D); SUNcGroup::Operations::DU(V5,V6,V5DV6);
            SUNcGroup::Operations::UU(V1DV2D,V3,V1DV2DV3); SUNcGroup::Operations::UU(V4,V5DV6,V4V5DV6);
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2DV3,V4V5DV6,V1DV2DV3V4V5DV6);
            
        }
        
        //V1D.V2.V3.V4D.V5D.V6
        void DUUDDU(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,SU_Nc_FUNDAMENTAL_FORMAT *V5,SU_Nc_FUNDAMENTAL_FORMAT *V6,SU_Nc_FUNDAMENTAL_FORMAT *V1DV2V3V4DV5DV6){
            
            //CREATE BUFFERS
            SU_Nc_FUNDAMENTAL_FORMAT V1DV2[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V5DV6[SUNcGroup::MatrixSize];SU_Nc_FUNDAMENTAL_FORMAT V1DV2V3[SUNcGroup::MatrixSize]; SU_Nc_FUNDAMENTAL_FORMAT V4DV5DV6[SUNcGroup::MatrixSize];
            
            //COMPUTE INDIVIDUAL PRODUCTS
            SUNcGroup::Operations::DU(V1,V2,V1DV2); SUNcGroup::Operations::DU(V5,V6,V5DV6);
            SUNcGroup::Operations::UU(V1DV2,V3,V1DV2V3); SUNcGroup::Operations::DU(V4,V5DV6,V4DV5DV6);
            
            
            //COMPUTE FINAL RESULTS
            SUNcGroup::Operations::UU(V1DV2V3,V4DV5DV6,V1DV2V3V4DV5DV6);
            
        }
        
        //V1+V2+V3+V4
        void ReTrITaSum(SU_Nc_FUNDAMENTAL_FORMAT *V1,SU_Nc_FUNDAMENTAL_FORMAT *V2,SU_Nc_FUNDAMENTAL_FORMAT *V3,SU_Nc_FUNDAMENTAL_FORMAT *V4,DOUBLE *SumProjection){
            
            //COMPUTE DIRECT SUM OF ALL ELEMENTS
            SU_Nc_FUNDAMENTAL_FORMAT MatrixSum[SUNcGroup::MatrixSize];
            
            for(int s=0;s<SUNcGroup::MatrixSize;s++){
                MatrixSum[s]=V1[s]+V2[s]+V3[s]+V4[s];
            }
            
            //COMPUTE TRACE PROJECTIONS
            SUNcGroup::Operations::ReTrIGenU(MatrixSum,SumProjection);
            
        }
        
        // COMPUTE EMid= -i SUNcLog(UNew.UOldDagger) AND UMid=exp( i EMid*c) UOld //
        void GeodesicInterpolation(DOUBLE c,SU_Nc_FUNDAMENTAL_FORMAT *UOld,SU_Nc_FUNDAMENTAL_FORMAT *UNew,SU_Nc_FUNDAMENTAL_FORMAT *UMid,SU_Nc_ALGEBRA_FORMAT *EMid){
            
            SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];
            
            SUNcGroup::Operations::UD(UNew,UOld,Buffer);
            SUNcAlgebra::Operations::MatrixILog(-DOUBLE(1.0),Buffer,EMid);
            
            SUNcAlgebra::Operations::MatrixIExp(-DOUBLE(c),EMid,Buffer);
            SUNcGroup::Operations::UU(Buffer,UOld,UMid);
            
        }
    
        ////////////////////////
        // BLOCKED MATRIX SUM //
        ////////////////////////
        /*
        void BlockRGTypeIMatrixSum(const SU_Nc_FUNDAMENTAL_FORMAT *VLine,const SU_Nc_FUNDAMENTAL_FORMAT *VStapleUp1,const SU_Nc_FUNDAMENTAL_FORMAT *VStapleDown1,const SU_Nc_FUNDAMENTAL_FORMAT *VStapleUp2,const SU_Nc_FUNDAMENTAL_FORMAT *VStapleDown2, SU_Nc_FUNDAMENTAL_FORMAT *Sum){
            
            // CREATE SUM //
            for(int s=0;s<SUNcGroup::MatrixSize;s++){
                //Sum[s]=(1.0-4.0/7.0)*VLine[s]+(1.0/7.0)*(VStapleUp1[s]+VStapleDown1[s]+VStapleUp2[s]+VStapleDown2[s]);
                Sum[s]=DOUBLE(1.0/3.0)*VLine[s]+(1.0/6.0)*(VStapleUp1[s]+VStapleDown1[s]+VStapleUp2[s]+VStapleDown2[s]);

            }
            
            // RE-UNITARIZE //
            SUNcGroup::Extended::MaxTraceProjection(Sum);
            
            
        }
         */
    }
    
}

namespace SUNcAlgebra {

    namespace Operations{
        
        DOUBLE ScalarProduct(SU_Nc_ALGEBRA_FORMAT *A1,SU_Nc_ALGEBRA_FORMAT *A2){
            
            DOUBLE A1DotA2=0.0;
            
            for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                A1DotA2+=A1[a]*A2[a];
            }
            
            return A1DotA2;
            
        }
    }
    
    namespace IO{
        
        std::string VectorToString(SU_Nc_ALGEBRA_FORMAT *alpha){
            
            std::stringstream sstm;
            
            sstm.precision(OUTPUT_PRECISION);
            
            for(int a=0;a<VectorSize;a++){
                sstm << alpha[a] << " ";
            }
            
            return sstm.str();
            
        }
    }
}

#endif
