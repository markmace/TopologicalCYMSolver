#ifndef _MINKOWSKI_WAVEVECTORS_CPP_
#define _MINKOWSKI_WAVEVECTORS_CPP_

namespace  PolarizationVectors {
    
    //BASIS VECTOR -- SHOULD NOT BE COLLINEAR TO ANY LATTICE MOMENTUM
    const DOUBLE BaseVector[3]={DOUBLE(1.0)/(D_SQRT3),DOUBLE(1.0)/(D_SQRT2),DOUBLE(1.0)/(D_SQRT2*D_SQRT3)};
    
    void Compute(INT pXIndex,COMPLEX cDpx,INT pYIndex,COMPLEX cDpy,INT pZIndex,COMPLEX cDpz,COMPLEX *Xi1,COMPLEX *Xi1Dot,COMPLEX *Xi2,COMPLEX *Xi2Dot){
        
        // MOMENTA AND FREQUENCY //
        DOUBLE pSqr=SQR_ABS(cDpx)+SQR_ABS(cDpy)+SQR_ABS(cDpz);    DOUBLE Omega=sqrt(pSqr);
        
        // SCALAR PRODUCTS //
        DOUBLE BaseVectorSqr=BaseVector[0]*BaseVector[0]+BaseVector[1]*BaseVector[1]+BaseVector[2]*BaseVector[2];
        
        COMPLEX BaseVectorDotP=BaseVector[0]*cDpx+BaseVector[1]*cDpy+BaseVector[2]*cDpz;
        
        DOUBLE Norm=sqrt(abs(BaseVectorSqr-BaseVectorDotP*conj(BaseVectorDotP)/pSqr));
        
        
        // SET POLARIZATION VECTORS //
        if(pSqr>DOUBLE(0.0)){
            
            //FIRST POLARIZATION VECTOR (v1=(e1 - (e1.p) p/ p^2)/sqrt(e1.e1-e1.p1/p^2) -- DIVIDED BY |P|
            Xi1[0]=(BaseVector[0]-BaseVectorDotP*conj(cDpx)/pSqr)/(Norm*sqrt(DOUBLE(2.0)*Omega));
            Xi1[1]=(BaseVector[1]-BaseVectorDotP*conj(cDpy)/pSqr)/(Norm*sqrt(DOUBLE(2.0)*Omega));
            Xi1[2]=(BaseVector[2]-BaseVectorDotP*conj(cDpz)/pSqr)/(Norm*sqrt(DOUBLE(2.0)*Omega));
            
            //SECOND POLARIZATION VECTOR (v2=p x e1)
            Xi2[0]=(cDpy*conj(Xi1[2])-cDpz*conj(Xi1[1]))/Omega;
            Xi2[1]=(cDpz*conj(Xi1[0])-cDpx*conj(Xi1[2]))/Omega;
            Xi2[2]=(cDpx*conj(Xi1[1])-cDpy*conj(Xi1[0]))/Omega;
                        
        }
        
        else{
            
            Xi1[0]=0.0;
            Xi1[1]=0.0;
            Xi1[2]=0.0;
            
            Xi2[0]=0.0;
            Xi2[1]=0.0;
            Xi2[2]=0.0;
            
        }
        
        // SET TIME DERIVATIVES OF POLARIZATION VECTORS //
        
        Xi1Dot[0]=ComplexI*Omega*Xi1[0];    Xi1Dot[1]=ComplexI*Omega*Xi1[1];   Xi1Dot[2]=ComplexI*Omega*Xi1[2];
        Xi2Dot[0]=ComplexI*Omega*Xi2[0];    Xi2Dot[1]=ComplexI*Omega*Xi2[1];   Xi2Dot[2]=ComplexI*Omega*Xi2[2];

        
        
    }
    
}
    
#endif
