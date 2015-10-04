#ifndef _SPECIAL_FUNCITONS_CPP_
#define _SPECIAL_FUNCITONS_CPP_

#include "SmartComplex.cpp"

/////////////////////////////////////////////////////////
// SPECIAL FUNCTION ROUTINES
/////////////////////////////////////////////////////////

namespace SpecialFunctions {
    
    DOUBLE Factorial (int n){
        
        if(n==0 || n==1){return 1;}
        else{return n*Factorial(n-1);}
    }
    
    
    //CONSTANTS FOR THE LANCZOS APPROXIMATION OF THE GAMMA FUNCTION FOR g=7 AND N=9 TRUNCATION OF THE SERIES
    static const int g=7;
    static const DOUBLE c[g+2]={0.99999999999980993, 676.5203681218851,-1259.1392167224028, 771.32342877765313, -176.61502916214059,12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,1.5056327351493116e-7};
    
    ///////////////////////////////////////////
    //GAMMA FUNCTION OF COMPLEX ARGUMENTS
    //////////////////////////////////////////
    
    SmartComplex gamma(COMPLEX z){
        
        //Use reflection formula if necessary
        //Gamma(1-z)Gamma(z)=pi/sin(pi z)
        
        if(real(z)<DOUBLE(0.5)){
            
            return PI/(SmartComplexFunctions::sine(PI*z)*gamma(DOUBLE(1.0)-z));
        }
        
        //Use Lanczos approximation to evaluate Gamma(z)
        //Gamma(z+1)=sqrt(2*pi)*(z+g+1/2)^(z+1/2)*exp(-(z+g+1/2)*A_g(z)
        
        //General expression for Ag
        //A_g(z)=0.5*p_0(g)+p_1(g)*(z/z+1)+p_2(g)*z*(z-1)/((z+1)*(z+2))+...
        
        //For fixed g this can be reexpressed as
        //A_g(z)=c_0+c_1/(z+1)+c_2/(z+2)+...
        
        COMPLEX Ag=c[0];
        for(int i=1;i<g+2;i++){
            Ag+=c[i]/(z+COMPLEX(-1+i,0));
        }
        
        //Define for convenience
        //t=(z-1+g+1/2)
        COMPLEX t=z+(g-DOUBLE(0.5));
        COMPLEX e=z-DOUBLE(0.5);
        
        //The final result is then given by
        //Gamma(z)=sqrt(2*pi)*exp(-t)*pow(t,e)*Ag
        
        //The complex power function contains the exp(pi Im(z)) factor
        //Re(x^a)=abs(x)^Re(a)*exp(-Im(a) arg(x))*cos(Im(a)log(abs(z))+arg(z)Re(a));
        //Im(x^a)=abs(x)^Re(a)*exp(-Im(a) arg(x))*sin(Im(a)log(abs(z))+arg(z)Re(a));
        
        //Store result in SmartComplex number
        SmartComplex GammaZ;
        
        DOUBLE resArg=imag(e)*log(abs(t))+arg(t)*real(e)+arg(Ag)-imag(t);
        
        GammaZ.number=sqrt(2.0*PI)*abs(Ag)*COMPLEX(cos(resArg),sin(resArg));
        GammaZ.order=real(e)*log(abs(t))-arg(t)*imag(e)-real(t);
        
        return GammaZ;
    }
    
    //////////////////////////////////////////////
    //BESSEL TYPE FUNTIONS OF COMPLEX ARGUMENTS
    //////////////////////////////////////////////
    
    
    SmartComplex BesselJSeries(COMPLEX nu,COMPLEX z){
        
        SmartComplex res;
        
        res.number=COMPLEX(0.0,0.0);
        res.order=0.0;
        
        for(int k=0;k<10;k++){
            res=res+pow(-1,k)*SmartComplexFunctions::power(DOUBLE(0.5)*z,COMPLEX(2*k,0)+nu)/(Factorial(k)*gamma(COMPLEX(k+1,0)+nu));
        }
        
        return res;
        
    }
    
    
    SmartComplex HankelH2Series(COMPLEX nu,COMPLEX z){
        
        //H2_a(x)=(J_-a(x) - exp(i pi a) J_a(x))/(-i sin(a pi))
        
        return (BesselJSeries(-nu,z)-SmartComplexFunctions::exponential(COMPLEX(0,1)*PI*nu)*BesselJSeries(nu,z))/(COMPLEX(0,-1)*SmartComplexFunctions::sine(PI*nu));
        
    }
    
    
    SmartComplex HankelH2DerivativeSeries(COMPLEX nu,COMPLEX z){
        
        //H2'_a(x)=(H2_a-1(x)-H2_a+1(x))/2;
        return DOUBLE(0.5)*(HankelH2Series(nu-DOUBLE(1.0),z)-HankelH2Series(nu+DOUBLE(1.0),z));
        
    }
    
}

#endif
