#ifndef __GAUGE_DEVAITION__CPP__
#define __GAUGE_DEVAITION__CPP__

#define SET_GAUGE_DEVIATION_BUFFERS() \
\
SU_Nc_FUNDAMENTAL_FORMAT w[SUNcGroup::MatrixSize];\
\
SU_Nc_FUNDAMENTAL_FORMAT U0G[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U1G[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U2G[SUNcGroup::MatrixSize];\
\
SU_Nc_FUNDAMENTAL_FORMAT U0DGxM[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U1DGyM[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT U2DGzM[SUNcGroup::MatrixSize];\
\
SU_Nc_ALGEBRA_FORMAT DeviationBuffer[SUNcAlgebra::VectorSize];

#define SET_GAUGE_UPDATE_BUFFERS() \
\
SU_Nc_FUNDAMENTAL_FORMAT wTilde[SUNcGroup::MatrixSize];\
SU_Nc_FUNDAMENTAL_FORMAT gNew[SUNcGroup::MatrixSize];

#define GET_LOCAL_GAUGE_VIOLATION(x,y,z) \
\
/*GET LOCAL GAUGE LINKS*/ \
COPY_SUNcMatrix(U0G,UGaugeTransformed->Get(x,y,z,0));\
COPY_SUNcMatrix(U1G,UGaugeTransformed->Get(x,y,z,1));\
COPY_SUNcMatrix(U2G,UGaugeTransformed->Get(x,y,z,2));\
\
/*GET INVERSE OF NEIGHBORING GAUGE LINKS*/\
SUNcGroup::Operations::Inverse(UGaugeTransformed->Get(x-1,y,z,0),U0DGxM);\
SUNcGroup::Operations::Inverse(UGaugeTransformed->Get(x,y-1,z,1),U1DGyM);\
SUNcGroup::Operations::Inverse(UGaugeTransformed->Get(x,y,z-1,2),U2DGzM);\
\
/*SET w(x)= g^{\mu\nu} U_{\nu}^{G}(x)+U^{(G)\dagger}(x+\nuhat)*/\
for(int alpha=0;alpha<SUNcGroup::MatrixSize;alpha++){\
    w[alpha]=Dynamics::gUpMetric[0]*SQR(Lattice::aScale/UGaugeTransformed->a[0])*(U0G[alpha]+U0DGxM[alpha])+Dynamics::gUpMetric[1]*SQR(Lattice::aScale/UGaugeTransformed->a[1])*(U1G[alpha]+U1DGyM[alpha])+(UGaugeTransformed->N[2]>1)*Dynamics::gUpMetric[2]*SQR(Lattice::aScale/UGaugeTransformed->a[2])*(U2G[alpha]+U2DGzM[alpha]);\
}\
\
/*COMPUTE COLOR TRACES TO MEASURE LOCAL DEVIATION*/\
SUNcGroup::Operations::ReTrIGenU(w,DeviationBuffer);

#define MONITOR_GAUGE_DEVIATION() \
    \
    /*MONITOR GAUGE DEVIATION*/\
    for(int a=0;a<SUNcAlgebra::VectorSize;a++){\
        MaxDeviation=std::max(MaxDeviation,DABS(DeviationBuffer[a])); \
        SqrSumDeviation+=SQR(DeviationBuffer[a]);\
    }\
    \
    /*MONITOR GAUGE FIXING FUNCTIONAL*/\
    DOUBLE LocalGaugeFunctional=SUNcGroup::Operations::ReTr(w); \
    MinGaugeFunctional=std::min(MinGaugeFunctional,LocalGaugeFunctional);\
    SumGaugeFunctional+=LocalGaugeFunctional;


#define MONITOR_GAUGE_UPDATE() \
    \
    /*CHECK UNITARITY*/\
    if(SUNcGroup::Operations::UnitarityNorm(wTilde)>std::pow(10.0,-MAX_DIGITS_PRECISION+2)){\
    \
    std::cerr << "################################################" << std::endl;\
    std::cerr << "##### UNITARITY VIOLATION -- GAUGE UPDATE ######" << std::endl;\
    std::cerr << SUNcGroup::Operations::UnitarityNorm(wTilde) << std::endl;\
    std::cerr << "################################################" << std::endl;\
    \
    exit(0);\
    }\
    \
    /*MONITOR STEP SIZE*/\
    DOUBLE LocalStepSize=abs(SUNcGroup::Operations::trIDMinusU(wTilde));\
    MaxStepSize=std::max(MaxStepSize,LocalStepSize); \
    SqrSumStepSize+=SQR(LocalStepSize);


#endif
