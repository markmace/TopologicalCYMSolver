#ifndef __CALC__PLAQUETTES__CPP__
#define __CALC__PLAQUETTES__CPP__


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                    ELEMENTARY PLAQUETTES                                                //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                 ------<------------<--------   ^                                        //
//                                                 |            ||            |  y|                                        //
//                                                 |            ||            |   |                                        //
//                                                 v U_{-xy}(x) ^v  U_{xy}(x) ^   |                                        //
//                                                 |            ||            |   ------->                                 //
//                                                 |            ||            |         x                                  //
//                                                 |            ||            |                                            //
//                                                 ------>------OO----->-------                                            //
//                                                 ------<------OO----<--------                                            //
//                                                 |            ||            |                                            //
//                                                 |            ||            |                                            //
//                                                 |U_{-x-y}(x) ^| U_{x-y}(x) ^                                            //
//                                                 v            |v            |                                            //
//                                                 |            ||            |                                            //
//                                                 |            ||            |                                            //
//                                                 ------>-------------->------                                            //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                            CONVENTIONS                                                       ///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                         //
//          U_{\mu\nu}(x)=U_{\mu}(x)U_{\nu}(x+\muhat)U_{\mu}^{\dagger}(x+\nuhat)U_{\nu}^{\dagger}(x)                       //
//         U_{-\mu\nu}(x)=U_{\nu}(x)U_{\mu}^{\dagger}(x+\nuhat-\muhat)U_{\nu}^{\dagger}(x-\muhat)U_{\mu}(x-\muhat)         //
//         U_{\mu-\nu}(x)=U_{\nu}^{\dagger}(x-\nuhat)U_{\mu}(x-\nuhat)U_{\nu}(x-\nuhat+\muhat)U_{\mu}^{\dagger}(x)         //
//         U_{-\mu-\nu}(x)=U_{\mu}^{\dagger}(x-\muhat)U_{\nu}^{\dagger}(x-\muhat)U_{\mu}(x-\muhat-\nuhat)U_{\nu}(x-\nuhat) //
//                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//BUFFERS FOR THE COMPUTATION OF ELEMENTARY PLAQUETTES
#define SET_ELEMENTARY_PLAQUETTE_BUFFERS() \
\
SU_Nc_FUNDAMENTAL_FORMAT Uzx[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT Uxy[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT Uyz[SUNcGroup::MatrixSize];

#define SET_NEIGHBORING_PLAQUETTE_BUFFERS() \
\
SU_Nc_FUNDAMENTAL_FORMAT UMzx[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UMxy[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UMyz[SUNcGroup::MatrixSize]; \
\
SU_Nc_FUNDAMENTAL_FORMAT UzMx[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UyMz[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UxMy[SUNcGroup::MatrixSize];

#define SET_EXTENDED_PLAQUETTE_BUFFERS() \
SU_Nc_FUNDAMENTAL_FORMAT UMyMz[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UMxMy[SUNcGroup::MatrixSize]; \
SU_Nc_FUNDAMENTAL_FORMAT UMzMx[SUNcGroup::MatrixSize];


//BUFFERS FOR THE COMPUTATION OF THE REAL PART OF THE TRACE iT^{a} TIMES THE ELEMENTARY PLAQUETTES
#define SET_ELEMENTARY_COLOR_TRACE_BUFFERS() \
\
DOUBLE ReTrITaUxy[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUzx[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUyz[SUNcAlgebra::VectorSize];

#define SET_NEIGHBORING_COLOR_TRACE_BUFFERS() \
\
DOUBLE ReTrITaUMxy[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUMyz[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUMzx[SUNcAlgebra::VectorSize]; \
\
DOUBLE ReTrITaUxMy[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUyMz[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUzMx[SUNcAlgebra::VectorSize];

#define SET_EXTENDED_COLOR_TRACE_BUFFERS() \
\
DOUBLE ReTrITaUMxMy[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUMyMz[SUNcAlgebra::VectorSize]; \
DOUBLE ReTrITaUMzMx[SUNcAlgebra::VectorSize];

//COMPUTATION OF ELEMENTARY PLAQUETTES
#define COMPUTE_ELEMENTARY_PLAQUETTES(x,y,z) \
\
SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,0),U->Get(x+1,y,z,1),U->Get(x,y+1,z,0),U->Get(x,y,z,1),Uxy); \
SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,1),U->Get(x,y+1,z,2),U->Get(x,y,z+1,1),U->Get(x,y,z,2),Uyz); \
SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,2),U->Get(x,y,z+1,0),U->Get(x+1,y,z,2),U->Get(x,y,z,0),Uzx); \

#define COMPUTE_NEIGHBORING_PLAQUETTES(x,y,z) \
\
SUNcGroup::AdvancedOperations::UDDU(U->Get(x,y,z,1),U->Get(x-1,y+1,z,0),U->Get(x-1,y,z,1),U->Get(x-1,y,z,0),UMxy); \
SUNcGroup::AdvancedOperations::UDDU(U->Get(x,y,z,2),U->Get(x,y-1,z+1,1),U->Get(x,y-1,z,2),U->Get(x,y-1,z,1),UMyz); \
SUNcGroup::AdvancedOperations::UDDU(U->Get(x,y,z,0),U->Get(x+1,y,z-1,2),U->Get(x,y,z-1,0),U->Get(x,y,z-1,2),UMzx); \
\
SUNcGroup::AdvancedOperations::DUUD(U->Get(x,y-1,z,1),U->Get(x,y-1,z,0),U->Get(x+1,y-1,z,1),U->Get(x,y,z,0),UxMy); \
SUNcGroup::AdvancedOperations::DUUD(U->Get(x,y,z-1,2),U->Get(x,y,z-1,1),U->Get(x,y+1,z-1,2),U->Get(x,y,z,1),UyMz); \
SUNcGroup::AdvancedOperations::DUUD(U->Get(x-1,y,z,0),U->Get(x-1,y,z,2),U->Get(x-1,y,z+1,0),U->Get(x,y,z,2),UzMx); \

#define COMPUTE_EXTENDED_PLAQUETTES(x,y,z) \
\
SUNcGroup::AdvancedOperations::DDUU(U->Get(x-1,y,z,0),U->Get(x-1,y-1,z,1),U->Get(x-1,y-1,z,0),U->Get(x,y-1,z,1),UMxMy); \
SUNcGroup::AdvancedOperations::DDUU(U->Get(x,y-1,z,1),U->Get(x,y-1,z-1,2),U->Get(x,y-1,z-1,1),U->Get(x,y,z-1,2),UMyMz); \
SUNcGroup::AdvancedOperations::DDUU(U->Get(x,y,z-1,2),U->Get(x-1,y,z-1,0),U->Get(x-1,y,z-1,2),U->Get(x-1,y,z,0),UMzMx);


//COMPUTATION OF THE PLAQUETTE TRACES
#define COMPUTE_ELEMENTARY_COLOR_TRACES() \
\
SUNcGroup::Operations::ReTrIGenU(Uxy,ReTrITaUxy); \
SUNcGroup::Operations::ReTrIGenU(Uyz,ReTrITaUyz); \
SUNcGroup::Operations::ReTrIGenU(Uzx,ReTrITaUzx);


#define COMPUTE_NEIGHBORING_COLOR_TRACES() \
\
SUNcGroup::Operations::ReTrIGenU(UMxy,ReTrITaUMxy); \
SUNcGroup::Operations::ReTrIGenU(UMyz,ReTrITaUMyz); \
SUNcGroup::Operations::ReTrIGenU(UMzx,ReTrITaUMzx); \
\
SUNcGroup::Operations::ReTrIGenU(UxMy,ReTrITaUxMy); \
SUNcGroup::Operations::ReTrIGenU(UyMz,ReTrITaUyMz); \
SUNcGroup::Operations::ReTrIGenU(UzMx,ReTrITaUzMx);

#define COMPUTE_EXTENDED_COLOR_TRACES() \
\
SUNcGroup::Operations::ReTrIGenU(UMxMy,ReTrITaUMxMy); \
SUNcGroup::Operations::ReTrIGenU(UMyMz,ReTrITaUMyMz); \
SUNcGroup::Operations::ReTrIGenU(UMzMx,ReTrITaUMzMx);




#endif
