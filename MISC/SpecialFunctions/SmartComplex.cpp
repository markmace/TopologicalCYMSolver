#ifndef _SMARTCOMPLEXCPP_
#define _SMARTCOMPLEXCPP_

////////////////////////////////////////////
// DATA FORMAT FOR SMARTCOMPLEX NUMBERS
////////////////////////////////////////////

class SmartComplex {
    
	
	//Smart Complex number saves a numer s in the format
	//s=exp(order)*val
	
public:
	
	//Order one complex number
	COMPLEX number;
	
	//Exponent of the number
	DOUBLE order;
	
	
	COMPLEX value(){
		return number*exp(order);
	}
	
    std::string toString(){
		
		DOUBLE FullExponent=order/log(DOUBLE(10.0));
		INT IntExponent=floor(FullExponent);
		
		COMPLEX val=number*pow(DOUBLE(10.0),FullExponent-IntExponent);
		
		
        std::stringstream ss;
        ss << "(" << real(val) << "," << imag(val) << ")E" << IntExponent;
		
		return ss.str();
	}
	
	//Constructor
	SmartComplex(){};
	
	//Destructor
	~SmartComplex(){};
	
};


////////////////////////////////////////////////
//OPERATIONS ON SMARTCOMPLEX NUMBERS
////////////////////////////////////////////////

//MULTIPLICATION
SmartComplex operator*(const SmartComplex  &num1, const SmartComplex &num2){
	
	SmartComplex res;
	
	COMPLEX num=num1.number*num2.number;
	INT expCorr=floor(log(abs(num)));
	
	res.number=num*exp(-expCorr);
	res.order=num1.order+num2.order+expCorr;
	
	return res;
	
}

//DIVISION
SmartComplex operator/(const SmartComplex  &num1, const SmartComplex &num2){
	
	SmartComplex res;
	
	COMPLEX num=num1.number/num2.number;
	INT expCorr=floor(log(abs(num)));
	
	res.number=num*exp(-expCorr);
	res.order=num1.order-num2.order+expCorr;
	
	return res;
	
}

//ADDITION
SmartComplex operator+(const SmartComplex  &num1, const SmartComplex &num2){
	
	SmartComplex res;
	
	if(num1.order>num2.order){
		
		COMPLEX num=num1.number+num2.number*exp(num2.order-num1.order);
		
		INT expCorr=floor(log(abs(num)));
		
		res.number=num*exp(-expCorr);
		res.order=num1.order+expCorr;
		
	}
	
	else{
		
		COMPLEX num=num2.number+num1.number*exp(num1.order-num2.order);
		
		int expCorr=floor(log(abs(num)));
		
		res.number=num*exp(-expCorr);
		res.order=num2.order+expCorr;
	}
	
	return res;
	
}

//SUBTRACTION
SmartComplex operator-(const SmartComplex  &num1, const SmartComplex &num2){
	
	SmartComplex res;
	
	if(num1.order>num2.order){
		
		COMPLEX num=num1.number-num2.number*exp(num2.order-num1.order);
		
		INT expCorr=floor(log(abs(num)));
		
		res.number=num*exp(-expCorr);
		res.order=num1.order+expCorr;
	}
	
	else{
		
		COMPLEX num=num1.number*exp(num1.order-num2.order)-num2.number;
		
		INT expCorr=floor(log(abs(num)));
		
		res.number=num*exp(-expCorr);
		res.order=num2.order+expCorr;
	}
	
	return res;
	
}

//SCALAR MULTIPLICATION WITH DOUBLE
SmartComplex operator*(const DOUBLE  &s, const SmartComplex &num){
	
	SmartComplex res;
	
	INT expCorr=floor(log(abs(s)));
	DOUBLE factor=s*exp(-expCorr);
	
	res.number=factor*num.number;
	res.order=expCorr+num.order;
	
	return res;
	
}

//SCALAR MULTIPLICATION WITH COMPLEX
SmartComplex operator*(const COMPLEX  &s, const SmartComplex &num){
	
	SmartComplex res;
	
	INT expCorr=floor(log(abs(s)));
	COMPLEX factor=s*exp(-expCorr);
	
	res.number=factor*num.number;
	res.order=expCorr+num.order;
	
	return res;
	
}

//SCALAR MULTIPLICATION OF INVERSE WITH DOUBLE
SmartComplex operator/(const DOUBLE  &s, const SmartComplex &num){
	
	SmartComplex res;
	
	INT expCorr=floor(log(abs(s)));
	DOUBLE factor=s*exp(-expCorr);
	
	res.number=factor/num.number;
	res.order=expCorr-num.order;
	
	return res;
	
}

//SCALAR DIVISION BY DOUBLE
SmartComplex operator/(const SmartComplex &num,const DOUBLE  &s){
	
	SmartComplex res;
	
	INT expCorr=floor(log(abs(s)));
	DOUBLE factor=s*exp(-expCorr);
	
	res.number=num.number/factor;
	res.order=num.order-expCorr;
	
	return res;
	
}

/////////////////////////////////////////////////////////
//ELEMENTARY FUNCTIONS OF SMARTCOMPLEX FUNCTIONS
/////////////////////////////////////////////////////////


namespace SmartComplexFunctions {
	
	//Computes sin(z) in smart complex format
	SmartComplex sine(COMPLEX z){
		
		SmartComplex res;
		
		res.number=COMPLEX(sin(real(z)),tanh(imag(z))*cos(real(z)));
		
		//USE COMPLETE log(cosh(x)) ONLY FOR x<10
		if(abs(imag(z))<DOUBLE(10.0)){
			res.order=log(cosh(imag(z)));
		}
		
		//ELSE USE log(cosh(x))~Abs[x]-log[2] (ABSOLUTE ERROR IS SMALLER THAN 1e-9 for |x|>10)
		else {
			res.order=abs(imag(z))-log(DOUBLE(2.0));
		}
		
		return res;
	}
	
	//Computes exp(z) in smart complex format
	SmartComplex exponential(COMPLEX z){
		
		SmartComplex res;
		
		res.number=COMPLEX(cos(imag(z)),sin(imag(z)));
		res.order=real(z);
		
		return res;
	}
	
	//Computes x^a in smart complex format
	SmartComplex power(COMPLEX x,COMPLEX a){
		
		SmartComplex res;
		
		DOUBLE resArg=imag(a)*log(abs(x))+arg(x)*real(a);
		
		res.number=COMPLEX(cos(resArg),sin(resArg));
		res.order=real(a)*log(abs(x))-arg(x)*imag(a);
		
		return res;
	}
	
}

#endif
