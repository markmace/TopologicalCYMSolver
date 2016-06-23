#ifndef __CABBIBO_MARINARI_PROJECTION_CPP__
#define __CABBIBO_MARINARI_PROJECTION_CPP__

#define GET_FIRST_CABBIBO_MARINARI_MATRIX() \
\
DOUBLE OneOverDetA1=DOUBLE(1.0)/sqrt(SQR_ABS(conj(VNew[MatrixIndex(0,0)])+VNew[MatrixIndex(1,1)])+SQR_ABS(conj(VNew[MatrixIndex(1,0)])-VNew[MatrixIndex(0,1)])); \
\
A[MatrixIndex(0,0)]=OneOverDetA1*(conj(VNew[MatrixIndex(0,0)])+VNew[MatrixIndex(1,1)]);\
A[MatrixIndex(0,1)]=OneOverDetA1*(conj(VNew[MatrixIndex(1,0)])-VNew[MatrixIndex(0,1)]);\
\
A[MatrixIndex(1,0)]=OneOverDetA1*(conj(VNew[MatrixIndex(0,1)])-VNew[MatrixIndex(1,0)]);\
A[MatrixIndex(1,1)]=OneOverDetA1*(conj(VNew[MatrixIndex(1,1)])+VNew[MatrixIndex(0,0)]);\
\
A[MatrixIndex(0,2)]=0.0;\
A[MatrixIndex(1,2)]=0.0;\
\
A[MatrixIndex(2,0)]=0.0;\
A[MatrixIndex(2,1)]=0.0;\
A[MatrixIndex(2,2)]=1.0;

#define GET_SECOND_CABBIBO_MARINARI_MATRIX() \
\
DOUBLE OneOverDetA2=DOUBLE(1.0)/sqrt(SQR_ABS(conj(VNew[MatrixIndex(0,0)])+VNew[MatrixIndex(2,2)])+SQR_ABS(conj(VNew[MatrixIndex(2,0)])-VNew[MatrixIndex(0,2)])); \
\
A[MatrixIndex(0,0)]=OneOverDetA2*(conj(VNew[MatrixIndex(0,0)])+VNew[MatrixIndex(2,2)]);\
A[MatrixIndex(0,2)]=OneOverDetA2*(conj(VNew[MatrixIndex(2,0)])-VNew[MatrixIndex(0,2)]);\
\
A[MatrixIndex(2,0)]=OneOverDetA2*(conj(VNew[MatrixIndex(0,2)])-VNew[MatrixIndex(2,0)]);\
A[MatrixIndex(2,2)]=OneOverDetA2*(conj(VNew[MatrixIndex(2,2)])+VNew[MatrixIndex(0,0)]);\
\
A[MatrixIndex(0,1)]=0.0;\
A[MatrixIndex(2,1)]=0.0;\
\
A[MatrixIndex(1,2)]=0.0;\
A[MatrixIndex(1,1)]=1.0;\
A[MatrixIndex(1,0)]=0.0;


#define GET_THIRD_CABBIBO_MARINARI_MATRIX() \
\
DOUBLE OneOverDetA3=DOUBLE(1.0)/sqrt(SQR_ABS(conj(VNew[MatrixIndex(1,1)])+VNew[MatrixIndex(2,2)])+SQR_ABS(conj(VNew[MatrixIndex(2,1)])-VNew[MatrixIndex(1,2)])); \
\
A[MatrixIndex(1,1)]=OneOverDetA3*(conj(VNew[MatrixIndex(1,1)])+VNew[MatrixIndex(2,2)]);\
A[MatrixIndex(1,2)]=OneOverDetA3*(conj(VNew[MatrixIndex(2,1)])-VNew[MatrixIndex(1,2)]);\
\
A[MatrixIndex(2,1)]=OneOverDetA3*(conj(VNew[MatrixIndex(1,2)])-VNew[MatrixIndex(2,1)]);\
A[MatrixIndex(2,2)]=OneOverDetA3*(conj(VNew[MatrixIndex(2,2)])+VNew[MatrixIndex(1,1)]);\
\
A[MatrixIndex(1,0)]=0.0;\
A[MatrixIndex(2,0)]=0.0;\
\
A[MatrixIndex(0,0)]=1.0;\
A[MatrixIndex(0,1)]=0.0;\
A[MatrixIndex(0,2)]=0.0;\
\

#define SET_CABBIBO_MARINARI_BUFFERS()\
SU_Nc_FUNDAMENTAL_FORMAT A[SUNcGroup::MatrixSize];\
\
SU_Nc_FUNDAMENTAL_FORMAT G[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT VNew[SUNcGroup::MatrixSize];\
\
SU_Nc_FUNDAMENTAL_FORMAT Buffer[SUNcGroup::MatrixSize];


#define CABBIBO_MARINARI_STEP()\
\
GET_FIRST_CABBIBO_MARINARI_MATRIX();\
\
SUNcGroup::Operations::UU(A,G,Buffer);\
COPY_SUNcMatrix(G,Buffer);\
\
SUNcGroup::Operations::UU(A,VNew,Buffer);\
COPY_SUNcMatrix(VNew,Buffer);\
\
GET_SECOND_CABBIBO_MARINARI_MATRIX();\
\
SUNcGroup::Operations::UU(A,G,Buffer);\
COPY_SUNcMatrix(G,Buffer);\
\
SUNcGroup::Operations::UU(A,VNew,Buffer);\
COPY_SUNcMatrix(VNew,Buffer);\
\
GET_THIRD_CABBIBO_MARINARI_MATRIX();\
\
SUNcGroup::Operations::UU(A,G,Buffer);\
COPY_SUNcMatrix(G,Buffer);\
\
SUNcGroup::Operations::UU(A,VNew,Buffer);\
COPY_SUNcMatrix(VNew,Buffer);



#endif
