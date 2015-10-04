#ifndef __INDEXING__CPP__
#define __INDEXING__CPP__

//CHECKERBOARD PARITY
#define CHECKERBOARD_EVEN_FLAG 0
#define CHECKERBOARD_ODD_FLAG 1
#define CHECKERBOARD_ALL_FLAG -1

#define CHECKERBOARD_PARITY(x,y,z) (((x)+(y)+(z))%2)

//GENERAL 3D INDEX
//#define Index3D(x,y,z)  (MOD((x),U->N[0])+U->N[0]*(MOD((y),U->N[1])+U->N[1]*MOD((z),U->N[2])))

//GENERAL 2D INDEX
//#define Index2D(x,y)  (MOD((x),U->N[0])+U->N[0]*(MOD((y),U->N[1])))

//GAUGE LINK INDEX -- MACRO TO RETURN THE POSITION OF FIRST ELEMENT OF FUNDAMENTAL GAUGE LINK
//#define GaugeLinkIndex(x,y,z,mu)    (SUNcGroup::MatrixSize*((mu)+Lattice::Dimension*Index3D((x),(y),(z))))

//GAUGE TRANSFORMATION INDEX -- MACRO TO RETURN THE POSITION OF FIRST ELEMENT OF GAUGE TRANSFORMATION
//#define GaugeTransformationIndex(x,y,z) (SUNcGroup::MatrixSize*(Index3D((x),(y),(z))))

//GAUGE LINK INDEX -- MACRO TO RETURN THE POSITION OF A CERTAIN COMPONENT OF COLORED VECTOR FIELD
//#define VectorFieldIndex(x,y,z,mu,a) ((a)+SUNcAlgebra::VectorSize*((mu)+Lattice::Dimension*Index3D((x),(y),(z))))

//GAUSS VIOLATION INDEX -- MACRO TO RETURN THE POSITION OF A CERTAIN COMPONENT OF COLORED FIELD
//#define GaussViolationIndex(x,y,z,a) ((a)+SUNcAlgebra::VectorSize*(Index3D((x),(y),(z))))

#endif
