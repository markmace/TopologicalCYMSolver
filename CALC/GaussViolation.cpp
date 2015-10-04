#ifndef __CALC__GAUSSLAW__VIOLATION__CPP__
#define __CALC__GAUSSLAW__VIOLATION__CPP__

#define SET_ADJOINT_BUFFERS() \
SU_Nc_ADJOINT_FORMAT U0Adjx[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U1Adjy[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U2Adjz[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
\
SU_Nc_ADJOINT_FORMAT U0AdjxM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U1AdjyM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U2AdjzM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];

#define SET_ASQR_BUFFERS() \
SU_Nc_ADJOINT_FORMAT U0AdjxP[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U1AdjyP[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U2AdjzP[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
\
SU_Nc_ADJOINT_FORMAT U0AdjxMM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U1AdjyMM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];\
SU_Nc_ADJOINT_FORMAT U2AdjzMM[SUNcAlgebra::VectorSize*SUNcAlgebra::VectorSize];

#define COMPUTE_NEIGHBORING_ADJOINTS(x,y,z) \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,0),U0Adjx); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,1),U1Adjy); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,2),U2Adjz); \
\
SUNcGroup::Operations::GetAdjoint(U->Get(x-1,y,z,0),U0AdjxM); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y-1,z,1),U1AdjyM); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z-1,2),U2AdjzM);

//THIS NEEDS TO BE CHECKED -- WRONG!!//
#define COMPUTE_ASQR_ADJOINTS(x,y,z) \
SUNcGroup::Operations::GetAdjoint(U->Get(x+1,y,z,0),U0AdjxP); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y+1,z,1),U1AdjyP); \
SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z+1,2),U2AdjzP); \
\
SU_Nc_FUNDAMENTAL_FORMAT U0Temp[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U1Temp[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U2Temp[SUNcGroup::MatrixSize];\
\
SUNcGroup::Operations::DD(U->Get(x-1,y,z,0),U->Get(x-2,y,z,0),U0Temp); \
SUNcGroup::Operations::DD(U->Get(x,y-1,z,1),U->Get(x,y-2,z,1),U1Temp); \
SUNcGroup::Operations::DD(U->Get(x,y,z-1,2),U->Get(x,y,z-2,2),U2Temp); \
\
SUNcGroup::Operations::GetAdjoint(U0Temp,U0AdjxMM); \
SUNcGroup::Operations::GetAdjoint(U1Temp,U1AdjyMM); \
SUNcGroup::Operations::GetAdjoint(U2Temp,U2AdjzMM);

//PROBABLY NEEDS TO BE CHANGED
#define COMPUTE_GAUSS_VIOLATION(LocalViolation,x,y,z,a) \
\
LocalViolation=DOUBLE(0.0); \
\
for(int b=0;b<SUNcAlgebra::VectorSize;b++){ \
LocalViolation+=(DELTA(a,b)*E->Get(x,y,z,0,b)[0]-U0AdjxM[SUNcGroup::AdjIndex(b,a)]*E->Get(x-1,y,z,0,b)[0])+(DELTA(a,b)*E->Get(x,y,z,1,b)[0]-U1AdjyM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y-1,z,1,b)[0])+(DELTA(a,b)*E->Get(x,y,z,2,b)[0]-U2AdjzM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y,z-1,2,b)[0]); \
} \
\
LocalViolation=LocalViolation*(Lattice::aScale*SQR(Lattice::aScale)/(Dynamics::MetricDeterminant*(E->a[0]*E->a[1]*E->a[2])));

#endif
