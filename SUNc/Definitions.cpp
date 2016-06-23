#ifndef __SUNc__DEFINITIONS_CPP__
#define __SUNc__DEFINITIONS_CPP__

//COPYING GAUGE LINKS 
#define COPY_SUNcMatrix(Destination,Origin)    std::memcpy((Destination),(Origin),(sizeof(SU_Nc_FUNDAMENTAL_FORMAT)*SUNcGroup::MatrixSize))

//DEFINITIONS OF TYPES FOR U(1) GAUGE GROUP
#if SU_Nc_FLAG==U1_FLAG

#define Nc 1

typedef DOUBLE SU_Nc_ALGEBRA_FORMAT;
typedef COMPLEX SU_Nc_FUNDAMENTAL_FORMAT;
typedef DOUBLE SU_Nc_ADJOINT_FORMAT;

#include "U1/GroupOperations.cpp"
#include "U1/AlgebraOperations.cpp"

#endif

//DEFINITIONS OF TYPES FOR SU(2) GAUGE GROUP
#if SU_Nc_FLAG==SU2_FLAG

#define Nc 2

typedef DOUBLE SU_Nc_ALGEBRA_FORMAT;
typedef DOUBLE SU_Nc_FUNDAMENTAL_FORMAT;
typedef DOUBLE SU_Nc_ADJOINT_FORMAT;

#include "SU2/GroupOperations.cpp"
#include "SU2/AlgebraOperations.cpp"

#endif

//DEFINITIONS OF TYPES FOR SU(3) GAUGE GROUP
#if SU_Nc_FLAG==SU3_FLAG

#define Nc 3

typedef DOUBLE SU_Nc_ALGEBRA_FORMAT;
typedef COMPLEX SU_Nc_FUNDAMENTAL_FORMAT;
typedef DOUBLE SU_Nc_ADJOINT_FORMAT;

#include "SU3/util/Diagonalization.cpp"
#include "SU3/util/ComplexEigensystem.cpp"
#include "SU3/GroupOperations.cpp"
#include "SU3/AlgebraOperations.cpp"

#endif

//INCLUDE COMPOSITE MATRIX PRODUCTS
#include "AdvancedOperations.cpp"


#endif
