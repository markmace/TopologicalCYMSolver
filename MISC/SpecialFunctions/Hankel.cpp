#ifndef _HANKEL_CPP_
#define _HANKEL_CPP_

#include "SpecialFunctions.cpp"

namespace SpecialFunctions {
    
    //STARTING VALUE FOR ODE
    COMPLEX HankelH2StartValue(COMPLEX nu,COMPLEX z){
        
        //Calculate Normalized value including exp(pi nu/2) factor
        SmartComplex Value=SmartComplexFunctions::exponential(DOUBLE(0.5)*PI*imag(nu))*HankelH2Series(nu,z);
        
        //Cast to complex number
        return Value.number*exp(Value.order);
        
    }
    
    //STARTING DERIVATIVE FOR ODE
    COMPLEX HankelH2DerivativeStartValue(COMPLEX nu,COMPLEX z){
        
        //Calculate Normalized value including exp(pi nu/2) factor
        SmartComplex Value=SmartComplexFunctions::exponential(DOUBLE(0.5)*PI*imag(nu))*HankelH2DerivativeSeries(nu,z);
        
        //Cast to complex number
        return Value.number*exp(Value.order);
        
    }
    
    
    //DATA FORMAT FOR HANKEL FUNCTION AND ITS DERIVATIVE
    struct HankelPair {
        COMPLEX Value;
        COMPLEX Derivative;
    };
    
    //CALCULATES HANKEL FUNCTION AND DERIVATIVE USING STARTING VALUES FROM SMALL T EXPANSION AND SOLVING BESSELS EQUATION
    void HankelH2Evolve(DOUBLE nu,DOUBLE xStart,HankelPair &H2Start,DOUBLE xGoal,HankelPair &H2Goal){
        
        //SOLVE BESSEL EQUATION
        INT Steps=0;
        DOUBLE x,dx;
        COMPLEX f,E;
        
        //Set starting value
        x=xStart;
        f=H2Start.Value;
        E=H2Start.Derivative*xStart;
        
        
        //LEAPFROG SOLVER FOR BESSEL EQUATION WITH ADAPTIVE STEPWIDTH
        //USE HIGHER ORDER SOLVER TO IMPROVE PRECISION OF THE FINAL RESULT
        if(x<xGoal){
            while(x<xGoal){
                
                //ADAPTIVE STEP SIZE
                if(abs(nu)>DOUBLE(1.0)){
                    dx=DOUBLE(0.001)*std::min(x/abs(nu),DOUBLE(1.0));
                }
                else{
                    dx=DOUBLE(0.001)*std::min(x,DOUBLE(1.0));
                }
                
                //Leap-frog solver
                f+=(dx/x)*E;
                
                E-=dx*x*(DOUBLE(1.0)+(nu*nu)/(x*x))*f;
                
                x+=dx;
                
                //Step counter
                Steps++;
            }
        }
        
        else if(x>xGoal){
            while(x>xGoal){
                
                //ADAPTIVE STEP SIZE
                if(abs(nu)>DOUBLE(1.0)){
                    dx=-DOUBLE(0.001)*std::min(x/abs(nu),DOUBLE(1.0));
                }
                else{
                    dx=-DOUBLE(0.001)*std::min(x,DOUBLE(1.0));
                }
                
                //Leap-frog solver
                f+=(dx/x)*E;
                
                E-=dx*x*(DOUBLE(1.0)+(nu*nu)/(x*x))*f;
                
                x+=dx;
                
                //Step counter
                Steps++;
            }
        }
        
        //Final step
        dx=xGoal-x;
        
        //Leap-frog solver
        f+=(dx/x)*E;
        
        E-=dx*(x+(nu*nu)/x)*f;
        
        x+=dx;	
        
        //Set Result
        H2Goal.Value=f;
        H2Goal.Derivative=E/x;
        
        
    }
    
    
    //CALCULATES HANKEL FUNCTION AND DERIVATIVE USING STARTING VALUES FROM SMALL T EXPANSION AND SOLVING BESSELS EQUATION
    void HankelH2(DOUBLE nu,DOUBLE xGoal,HankelPair &H2Goal){
        
        //Starting point
        DOUBLE xStart=DOUBLE(1.0);
        
        HankelPair H2Start;
        
        //Initialize at xStart using Series Expansion for small x 
        if(nu!=0){
            
            H2Start.Value=HankelH2StartValue(COMPLEX(0,nu),COMPLEX(xStart,0));
            H2Start.Derivative=HankelH2DerivativeStartValue(COMPLEX(0,nu),COMPLEX(xStart,0));
        }
        
        //For nu=0 use fixed values at x=1 
        else {
            
            //STARTING VALUES AT x=1
            xStart=DOUBLE(1.0);
            
            //FIXED VALUES FROM MATHEMATICA
            H2Start.Value=COMPLEX(0.76519769,-0.08825696);
            H2Start.Derivative=COMPLEX(-0.44005059,-0.78121282);
        }
        
        //EVOLVE TO XGOAL
        HankelH2Evolve(nu,xStart,H2Start,xGoal,H2Goal);
        
    }
    
}

#endif
