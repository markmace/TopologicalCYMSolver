#ifndef __AVG_FIELD_STRENGTH__CPP__
#define __AVG_FIELD_STRENGTH__CPP__

////////////////////////////////////////////////////////
// BUFFERS FOR COMPUTATION OF AVERAGE FIELD STRENGTHS //
////////////////////////////////////////////////////////

// MAGNETIC //
#define SET_AVG_MAGNETIC_FIELD_STRENGTH_BUFFERS() \
\
SET_ELEMENTARY_PLAQUETTE_BUFFERS(); \
SET_NEIGHBORING_PLAQUETTE_BUFFERS(); \
SET_EXTENDED_PLAQUETTE_BUFFERS(); \
\
SU_Nc_FUNDAMENTAL_FORMAT UzxAvg[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UxyAvg[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UyzAvg[SUNcGroup::MatrixSize]; \
\
DOUBLE B0Loc[SUNcAlgebra::VectorSize];\
DOUBLE B1Loc[SUNcAlgebra::VectorSize];\
DOUBLE B2Loc[SUNcAlgebra::VectorSize];\
\
DOUBLE B0SqrLoc,B1SqrLoc,B2SqrLoc;


// ELECTRIC //
#define SET_AVG_ELECTRIC_FIELD_STRENGTH_BUFFERS() \
\
SET_ADJOINT_BUFFERS(); \
\
DOUBLE E0Loc[SUNcAlgebra::VectorSize];\
DOUBLE E1Loc[SUNcAlgebra::VectorSize];\
DOUBLE E2Loc[SUNcAlgebra::VectorSize];\
\
DOUBLE E0SqrLoc,E1SqrLoc,E2SqrLoc;

// ALL //
#define SET_AVG_FIELD_STRENGTH_BUFFERS() \
\
SET_AVG_MAGNETIC_FIELD_STRENGTH_BUFFERS();\
SET_AVG_ELECTRIC_FIELD_STRENGTH_BUFFERS();


/////////////////////////////////////////////////////////////
// COMPUTATION OF AVERAGE FIELD STRENGTHS AT POINT (X,Y,Z) //
/////////////////////////////////////////////////////////////

// MAGNETIC //
#define COMPUTE_AVG_MAGNETIC_FIELD_STRENGTH(x,y,z) \
\
COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z); \
COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z); \
COMPUTE_EXTENDED_PLAQUETTES(x,y,z); \
\
SUNcGroup::AvgMatrix(Uzx,UMzx,UzMx,UMzMx,UzxAvg); \
SUNcGroup::AvgMatrix(Uxy,UMxy,UxMy,UMxMy,UxyAvg); \
SUNcGroup::AvgMatrix(Uyz,UMyz,UyMz,UMyMz,UyzAvg); \
\
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UyzAvg,B0Loc); \
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UzxAvg,B1Loc); \
SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UxyAvg,B2Loc); \
\
B0SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UyzAvg); \
B1SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UzxAvg); \
B2SqrLoc=DOUBLE(4.0)*SUNcGroup::Operations::ReTrIDMinusU(UxyAvg);


// ELECTRIC //
#define COMPUTE_AVG_ELECTRIC_FIELD_STRENGTH(x,y,z) \
\
for(int a=0;a<SUNcAlgebra::VectorSize;a++){\
    E0Loc[a]=0.0;   E1Loc[a]=0.0;   E2Loc[a]=0.0;\
}\
\
COMPUTE_NEIGHBORING_ADJOINTS(x,y,z); \
\
for(int a=0;a<SUNcAlgebra::VectorSize;a++){\
    for(int b=0;b<SUNcAlgebra::VectorSize;b++){ \
        E0Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,0,b)[0]+U0AdjxM[SUNcGroup::AdjIndex(b,a)]*E->Get(x-1,y,z,0,b)[0]);\
        E1Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,1,b)[0]+U1AdjyM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y-1,z,1,b)[0]);\
        E2Loc[a]+=DOUBLE(0.5)*(DELTA(a,b)*E->Get(x,y,z,2,b)[0]+U2AdjzM[SUNcGroup::AdjIndex(b,a)]*E->Get(x,y,z-1,2,b)[0]);\
    } \
}\
\
E0SqrLoc=0.0;   E1SqrLoc=0.0;   E2SqrLoc=0.0; \
\
for(int a=0;a<SUNcAlgebra::VectorSize;a++){\
\
    E0SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,0,a)[0])+SQR(E->Get(x-1,y,z,0,a)[0])); \
    E1SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,1,a)[0])+SQR(E->Get(x,y-1,z,1,a)[0])); \
    E2SqrLoc+=DOUBLE(0.5)*(SQR(E->Get(x,y,z,2,a)[0])+SQR(E->Get(x,y,z-1,2,a)[0])); \
\
}

// ALL //
#define COMPUTE_AVG_FIELD_STRENGTH(x,y,z) \
\
COMPUTE_AVG_MAGNETIC_FIELD_STRENGTH(x,y,z);\
COMPUTE_AVG_ELECTRIC_FIELD_STRENGTH(x,y,z);


#endif
