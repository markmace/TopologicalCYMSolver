#ifndef __CALC__GAUSSLAW__VIOLATION__CPP__
#define __CALC__GAUSSLAW__VIOLATION__CPP__


// GAUSS LAW BUFFERS //
#define SET_GAUSS_LAW_BUFFERS() \
\
SU_Nc_ALGEBRA_FORMAT ExDown[SUNcAlgebra::VectorSize]; \
SU_Nc_ALGEBRA_FORMAT EyDown[SUNcAlgebra::VectorSize]; \
SU_Nc_ALGEBRA_FORMAT EzDown[SUNcAlgebra::VectorSize]; \

// GAUSS LAW VIOLATION //
#define COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z) \
\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x-1,y,z,0),E->Get(x-1,y,z,0,0),ExDown);\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y-1,z,1),E->Get(x,y-1,z,1,0),EyDown);\
SUNcAlgebra::Operations::InverseAdjointMultiplication(U->Get(x,y,z-1,2),E->Get(x,y,z-1,2,0),EzDown);\
\
DOUBLE GaussLawNormalizationFactor=(Lattice::aScale*SQR(Lattice::aScale)/(Dynamics::MetricDeterminant*(E->a[0]*E->a[1]*E->a[2])));\
\
for(INT a=0;a<SUNcAlgebra::VectorSize;a++){\
    LocalViolation[a]=GaussLawNormalizationFactor*((E->Get(x,y,z,0,a)[0]-ExDown[a]) + (E->Get(x,y,z,1,a)[0]-EyDown[a]) + (E->Get(x,y,z,2,a)[0]-EzDown[a]));\
}

#endif
