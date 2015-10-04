#ifndef _BJORKEN_WAVEVECTORS_CPP_
#define _BJORKEN_WAVEVECTORS_CPP_

#include "../MISC/SpecialFunctions/Hankel.cpp"

namespace PolarizationVectors {
    
    void Compute(INT pXIndex,COMPLEX cDpx,INT pYIndex,COMPLEX cDpy,INT nuIndex,COMPLEX cDNu,COMPLEX *Xi1,COMPLEX *Xi1Dot,COMPLEX *Xi2,COMPLEX *Xi2Dot){
        
        // DETERMINE MAGNITUDE OF MOMENTA //
        DOUBLE tau0=Dynamics::tau; DOUBLE nu=abs(cDNu);    DOUBLE pT=sqrt(SQR_ABS(cDpx)+SQR_ABS(cDpy));
        
        //ZERO TRANSVERSE MOMENTUM MODES
        if(pXIndex==0 && pYIndex==0){
            
            //ZERO MODE SET TO ZERO
            if(nuIndex==0 || nuIndex==U->N[2]/2){
                
                //FIRST SOLUTION
                Xi1[0]=0.0;
                Xi1[1]=0.0;
                Xi1[2]=0.0;
                
                Xi1Dot[0]=0.0;
                Xi1Dot[1]=0.0;
                Xi1Dot[2]=0.0;
                
                
                //SECOND SOLUTION
                Xi2[0]=0.0;
                Xi2[1]=0.0;
                Xi2[2]=0.0;
                
                Xi2Dot[0]=0.0;
                Xi2Dot[1]=0.0;
                Xi2Dot[2]=0.0;
                
            }
            
            
            //NON-ZERO NU MODES
            else {
                
                //Initial values
                COMPLEX xi0=COMPLEX(DOUBLE(1.0)/sqrt(DOUBLE(2.0)*abs(nu)),0);
                COMPLEX xiPrime0=COMPLEX(0,-abs(nu)/tau0)*xi0;
                
                //FIRST SOLUTION
                Xi1[0]=xi0;
                Xi1[1]=0.0;
                Xi1[2]=0.0;
                
                Xi1Dot[0]=xiPrime0;
                Xi1Dot[1]=0.0;
                Xi1Dot[2]=0.0;
                
                
                //SECOND SOLUTION
                Xi2[0]=0.0;
                Xi2[1]=xi0;
                Xi2[2]=0.0;
                
                Xi2Dot[0]=0.0;
                Xi2Dot[1]=xiPrime0;
                Xi2Dot[2]=0.0;
                
            }
        }
        
        //GENERIC MODES
        else {
            
            //Get Hankel function
            SpecialFunctions::HankelPair H2;
            
            SpecialFunctions::HankelH2(nu,pT*tau0,H2);
            
            //FIRST SOLUTION
            
            //Norm Factors for Xi and XiDot
            DOUBLE NormXi=DOUBLE(0.5)*sqrt(PI)/pT;
            DOUBLE NormXiDot=DOUBLE(0.5)*sqrt(PI);
            
            //Polarization Vector
            Xi1[0] = -NormXi*cDpy*H2.Value;
            Xi1[1] =  NormXi*cDpx*H2.Value;
            Xi1[2] =  COMPLEX(0.0,0.0);
            
            //Derivative of Polarization Vector
            Xi1Dot[0]= -NormXiDot*cDpy*H2.Derivative;
            Xi1Dot[1]=  NormXiDot*cDpx*H2.Derivative;
            Xi1Dot[2]=  COMPLEX(0.0,0.0);
            
            //SECOND SOLUTION
            
            //Norm Factor
            DOUBLE c2=DOUBLE(0.5)*tau0*sqrt(PI)*(pT*tau0);
            
            //R and R' at t=tau0
            COMPLEX Rtau0=-(pT*tau0)/(nu*nu+(pT*tau0)*(pT*tau0))*c2*H2.Derivative;
            COMPLEX RN0=Rtau0;
            
            COMPLEX RTPrime0=c2/tau0*H2.Value;
            COMPLEX RNPrime0=RTPrime0;
            
            //Polarization Vector
            Xi2[0] =  cDNu*conj(cDpx)/((pT*tau0)*(pT*tau0))*Rtau0;
            Xi2[1] =  cDNu*conj(cDpy)/((pT*tau0)*(pT*tau0))*Rtau0;
            Xi2[2] =  -RN0;
            
            //Derivative of Polarization Vector
            Xi2Dot[0]= cDNu*conj(cDpx)/((pT*tau0)*(pT*tau0))*RTPrime0;
            Xi2Dot[1]= cDNu*conj(cDpy)/((pT*tau0)*(pT*tau0))*RTPrime0;
            Xi2Dot[2]= -RNPrime0;
            
        }
		
    }
    
}

#endif
