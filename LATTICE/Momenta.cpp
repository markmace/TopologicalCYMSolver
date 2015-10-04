#ifndef __MOMENTA__CPP__
#define __MOMENTA__CPP__

#define GetMomenta(pXIndex,pYIndex,pZIndex,cDpx,cDpy,cDpz) \
\
cDpx=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pXIndex)/DOUBLE(GLinks::U->N[0])))/GLinks::U->a[0]; \
cDpy=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pYIndex)/DOUBLE(GLinks::U->N[1])))/GLinks::U->a[1]; \
cDpz=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pZIndex)/DOUBLE(GLinks::U->N[2])))/GLinks::U->a[2];

#define GetAbsMomentum(pXIndex,pYIndex,pZIndex,pAbs)\
\
cDpx=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pXIndex)/DOUBLE(GLinks::U->N[0])))/GLinks::U->a[0]; \
cDpy=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pYIndex)/DOUBLE(GLinks::U->N[1])))/GLinks::U->a[1]; \
cDpz=-ComplexI*(DOUBLE(1.0)-exp(-ComplexI*DOUBLE(2.0)*PI*DOUBLE(pZIndex)/DOUBLE(GLinks::U->N[2])))/GLinks::U->a[2]; \
\
pAbs=sqrt(SQR_ABS(cDpx)+SQR_ABS(cDpy)+SQR_ABS(cDpz)); 

#endif

