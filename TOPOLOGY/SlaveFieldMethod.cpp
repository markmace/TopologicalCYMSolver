namespace ChernSimonsNumber{
    
    namespace SlaveFieldMethod{
        
        // EACH POINT SHIFT IS A THREE DIMENSIONAL VECTOR //
        struct PointShift{
            
            INT X[Lattice::Dimension];
            
        };
        
        
        // EACH TETRAHEDRON CONSISTS OF FOUR POINT SHIFTS //
        struct Tetrahedron{
            
            PointShift P[4];
            
        };
        
        
        // EACH CELL CONSISTS OF FIVE TETRAHEDRA -- DIFFERENT FOR EVEN AND ODD CHECKERBOARD CELLS //
        namespace EvenCellDecomposition{
            
            Tetrahedron T[5];
            
            ///////////////////////////////////////////////////////////
            // DEFINITION OF TETRAHEADRA ACCORDING TO hep-ph/9703266 //
            ///////////////////////////////////////////////////////////
            
            void Set(){
                // DEFINITION OF THE FIRST TETRAHEDRON //
                T[0].P[0].X[0]=0;   T[0].P[0].X[1]=0;   T[0].P[0].X[2]=0;
                T[0].P[1].X[0]=1;   T[0].P[1].X[1]=0;   T[0].P[1].X[2]=0;
                T[0].P[2].X[0]=0;   T[0].P[2].X[1]=1;   T[0].P[2].X[2]=0;
                T[0].P[3].X[0]=0;   T[0].P[3].X[1]=0;   T[0].P[3].X[2]=1;
                
                // DEFINITION OF THE SECOND TETRAHEDRON //
                T[1].P[0].X[0]=1;   T[1].P[0].X[1]=0;   T[1].P[0].X[2]=0;
                T[1].P[1].X[0]=0;   T[1].P[1].X[1]=0;   T[1].P[1].X[2]=1;
                T[1].P[2].X[0]=1;   T[1].P[2].X[1]=0;   T[1].P[2].X[2]=1;
                T[1].P[3].X[0]=1;   T[1].P[3].X[1]=1;   T[1].P[3].X[2]=1;
                
                // DEFINITION OF THE THIRD TETRAHEDRON //
                T[2].P[0].X[0]=0;   T[2].P[0].X[1]=1;   T[2].P[0].X[2]=0;
                T[2].P[1].X[0]=1;   T[2].P[1].X[1]=0;   T[2].P[1].X[2]=0;
                T[2].P[2].X[0]=1;   T[2].P[2].X[1]=1;   T[2].P[2].X[2]=0;
                T[2].P[3].X[0]=1;   T[2].P[3].X[1]=1;   T[2].P[3].X[2]=1;
                
                // DEFINITION OF THE FOURTH TETRAHEDRON //
                T[3].P[0].X[0]=0;   T[3].P[0].X[1]=0;   T[3].P[0].X[2]=1;
                T[3].P[1].X[0]=0;   T[3].P[1].X[1]=1;   T[3].P[1].X[2]=0;
                T[3].P[2].X[0]=0;   T[3].P[2].X[1]=1;   T[3].P[2].X[2]=1;
                T[3].P[3].X[0]=1;   T[3].P[3].X[1]=1;   T[3].P[3].X[2]=1;
                
                // DEFINITION OF THE FIFTH TETRAHEDRON //
                T[4].P[0].X[0]=1;   T[4].P[0].X[1]=0;   T[4].P[0].X[2]=0;
                T[4].P[1].X[0]=0;   T[4].P[1].X[1]=1;   T[4].P[1].X[2]=0;
                T[4].P[2].X[0]=0;   T[4].P[2].X[1]=0;   T[4].P[2].X[2]=1;
                T[4].P[3].X[0]=1;   T[4].P[3].X[1]=1;   T[4].P[3].X[2]=1;
                
            }
            
        }
        
        namespace OddCellDecomposition {
            
            Tetrahedron T[5];
            
            void Set(){
                
                ///////////////////////////////////////////////////////////
                // DEFINITION OF TETRAHEADRA ACCORDING TO hep-ph/9703266 //
                ///////////////////////////////////////////////////////////
                
                // DEFINITION OF THE FIRST TETRAHEDRON //
                T[0].P[0].X[0]=1;   T[0].P[0].X[1]=0;   T[0].P[0].X[2]=0;
                T[0].P[1].X[0]=1;   T[0].P[1].X[1]=1;   T[0].P[1].X[2]=0;
                T[0].P[2].X[0]=0;   T[0].P[2].X[1]=0;   T[0].P[2].X[2]=0;
                T[0].P[3].X[0]=1;   T[0].P[3].X[1]=0;   T[0].P[3].X[2]=1;
                
                // DEFINITION OF THE SECOND TETRAHEDRON //
                T[1].P[0].X[0]=0;   T[1].P[0].X[1]=1;   T[1].P[0].X[2]=0;
                T[1].P[1].X[0]=1;   T[1].P[1].X[1]=1;   T[1].P[1].X[2]=0;
                T[1].P[2].X[0]=0;   T[1].P[2].X[1]=1;   T[1].P[2].X[2]=1;
                T[1].P[3].X[0]=0;   T[1].P[3].X[1]=0;   T[1].P[3].X[2]=0;
                
                // DEFINITION OF THE THIRD TETRAHEDRON //
                T[2].P[0].X[0]=0;   T[2].P[0].X[1]=0;   T[2].P[0].X[2]=1;
                T[2].P[1].X[0]=1;   T[2].P[1].X[1]=0;   T[2].P[1].X[2]=1;
                T[2].P[2].X[0]=0;   T[2].P[2].X[1]=0;   T[2].P[2].X[2]=0;
                T[2].P[3].X[0]=0;   T[2].P[3].X[1]=1;   T[2].P[3].X[2]=1;
                
                // DEFINITION OF THE FOURTH TETRAHEDRON //
                T[3].P[0].X[0]=1;   T[3].P[0].X[1]=1;   T[3].P[0].X[2]=1;
                T[3].P[1].X[0]=1;   T[3].P[1].X[1]=1;   T[3].P[1].X[2]=0;
                T[3].P[2].X[0]=1;   T[3].P[2].X[1]=0;   T[3].P[2].X[2]=1;
                T[3].P[3].X[0]=0;   T[3].P[3].X[1]=1;   T[3].P[3].X[2]=1;
                
                // DEFINITION OF THE FIFTH TETRAHEDRON //
                T[4].P[0].X[0]=0;   T[4].P[0].X[1]=0;   T[4].P[0].X[2]=0;
                T[4].P[1].X[0]=1;   T[4].P[1].X[1]=1;   T[4].P[1].X[2]=0;
                T[4].P[2].X[0]=0;   T[4].P[2].X[1]=1;   T[4].P[2].X[2]=1;
                T[4].P[3].X[0]=1;   T[4].P[3].X[1]=0;   T[4].P[3].X[2]=1;
                
            }
            
        }
        
        
        // GET PROJECTION ON THREE SPHERE //
        void GetSphereProjection(SU_Nc_FUNDAMENTAL_FORMAT *Q,DOUBLE *q){
            
            SUNcGroup::Operations::ReTrIGenU(Q,&q[0]);     q[3]=DOUBLE(0.5)*SUNcGroup::Operations::ReTr(Q);
            
        }
        
        // GET REFERENCE POINT FOR COMPUTATION OF THE DEGREE //
        void GetReferencePoint(DOUBLE *q){
            
            SU_Nc_FUNDAMENTAL_FORMAT Q[SUNcGroup::MatrixSize];
            
            RandomNumberGenerator::SUNcMatrix(DOUBLE(2.0)*D_SQRT2*PI,Q);
            
            GetSphereProjection(Q,q);
            
        }
        
        void Init(){
            
            EvenCellDecomposition::Set();
            OddCellDecomposition::Set();
            
        }
        
        DOUBLE LeviCivitaContraction(DOUBLE *p1,DOUBLE *p2,DOUBLE *p3,DOUBLE *p4){
            
            return p1[3]*p2[2]*p3[1]*p4[0] - p1[2]*p2[3]*p3[1]*p4[0] - p1[3]*p2[1]*p3[2]*p4[0] + p1[1]*p2[3]*p3[2]*p4[0] + p1[2]*p2[1]*p3[3]*p4[0] - p1[1]*p2[2]*p3[3]*p4[0] - p1[3]*p2[2]*p3[0]*p4[1] + p1[2]*p2[3]*p3[0]*p4[1] + p1[3]*p2[0]*p3[2]*p4[1] - p1[0]*p2[3]*p3[2]*p4[1] - p1[2]*p2[0]*p3[3]*p4[1] + p1[0]*p2[2]*p3[3]*p4[1] + p1[3]*p2[1]*p3[0]*p4[2] - p1[1]*p2[3]*p3[0]*p4[2] - p1[3]*p2[0]*p3[1]*p4[2] + p1[0]*p2[3]*p3[1]*p4[2] + p1[1]*p2[0]*p3[3]*p4[2] - p1[0]*p2[1]*p3[3]*p4[2] - p1[2]*p2[1]*p3[0]*p4[3] + p1[1]*p2[2]*p3[0]*p4[3] + p1[2]*p2[0]*p3[1]*p4[3] - p1[0]*p2[2]*p3[1]*p4[3] - p1[1]*p2[0]*p3[2]*p4[3] + p1[0]*p2[1]*p3[2]*p4[3];
            
        }
        
        INT CheckSign(DOUBLE x){
            
            if(DABS(x)<std::pow(10.0,-14)){
                return 0;
            }
            else{
                return SIGN(x);
            }
            
        }
        
        INT Orientation(DOUBLE *p1,DOUBLE *p2,DOUBLE *p3,DOUBLE *p4){
            
            return CheckSign(LeviCivitaContraction(p1,p2,p3,p4));
            
        }
        
        
        DOUBLE CheckOrientation(DOUBLE p[4][4],DOUBLE *q){
            
            DOUBLE OX=Orientation(p[0],p[1],p[2],p[3]);
            
            DOUBLE O0=Orientation(  q ,p[1],p[2],p[3]);
            DOUBLE O1=Orientation(p[0],  q ,p[2],p[3]);
            DOUBLE O2=Orientation(p[0],p[1],  q ,p[3]);
            DOUBLE O3=Orientation(p[0],p[1],p[2], q  );
            
            if(OX==O0 && OX==O1 && OX==O2 && OX==O3){
                return OX;
            }
            
            else{
                return 0;
            }
            
            
        }
        
        
        DOUBLE NCS(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeTransformations *G){
            
            // SPHERE REFERENCE POINT //
            DOUBLE qSphere[4];   
            
            // SPHERE TETRAHEDRON POINTS //
            DOUBLE pSphere[4][4];
            
            // GET REFERENCE POINT //
            GetReferencePoint(qSphere);
            
            // DETERMINE CORNERS OF THE TETRAHEDRON //
            INT xPos[Lattice::Dimension]; INT pXPos[5][4][3];
            
            // CHERN SIMONS NUMBER //
            DOUBLE ChernSimonsNumber=0.0;
            
            // SUM OVER SPACE GRID AT A GIVEN TIMESLICE //
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                                                
                        // SET POSITION //
                        xPos[0]=x; xPos[1]=y; xPos[2]=z;
                        
                        if(CHECKERBOARD_PARITY(x,y,z)==CHECKERBOARD_EVEN_FLAG){
                            
                            // USE EVEN CELL DECOMPOSITION //
                            using namespace EvenCellDecomposition;
                            
                            // GET ALL POSITIONS OF ALL TETRAHEDRA CORNERS //
                            for(INT t=0;t<5;t++){
                                for(INT p=0;p<4;p++){
                                    for(INT x=0;x<Lattice::Dimension;x++){
                                        pXPos[t][p][x]=xPos[x]+T[t].P[p].X[x];
                                    }
                                }
                            }
                            
                            // COMPUTE DEGREE OF THE MAP //
                            for(INT t=0;t<5;t++){
                                
                                // GET ALL VALUES ON THE THREE SPHERE AT TETRAHEDRON CORNERS //
                                for(INT p=0;p<4;p++){
                                    GetSphereProjection(G->Get(pXPos[t][p][0],pXPos[t][p][1],pXPos[t][p][2]),&pSphere[p][0]);
                                }
                                
                                // DETERMINE WHETHER Q IS CONTAINED IN THE INVERSE IMAGE OF THE MAP AND IF SO DETERMINE THE ORIENTATION //
                                ChernSimonsNumber+=CheckOrientation(pSphere,qSphere);
                                
                            }
                            
                        }
                        
                        if(CHECKERBOARD_PARITY(x,y,z)==CHECKERBOARD_ODD_FLAG){
                            
                            // USE ODD CELL DECOMPOSITION //
                            using namespace OddCellDecomposition;
                            
                            // GET ALL POSITIONS OF ALL TETRAHEDRA CORNERS //
                            for(INT t=0;t<5;t++){
                                for(INT p=0;p<4;p++){
                                    for(INT x=0;x<Lattice::Dimension;x++){
                                        pXPos[t][p][x]=xPos[x]+T[t].P[p].X[x];
                                    }
                                }
                            }
                            
                            // COMPUTE DEGREE OF THE MAP //
                            for(INT t=0;t<5;t++){
                                
                                // GET ALL VALUES ON THE THREE SPHERE AT TETRAHEDRON CORNERS //
                                for(INT p=0;p<4;p++){
                                    GetSphereProjection(G->Get(pXPos[t][p][0],pXPos[t][p][1],pXPos[t][p][2]),&pSphere[p][0]);
                                }
                                
                                // DETERMINE WHETHER Q IS CONTAINED IN THE INVERSE IMAGE OF THE MAP AND IF SO DETERMINE THE ORIENTATION //
                                ChernSimonsNumber+=CheckOrientation(pSphere,qSphere);
                                
                            }
                            
                        }
                        
                        
                    }
                }
            }
            
            //mstd::cerr << "#MEASURING SLAVE FIELD AT T=" << Dynamics::Time() << std::endl;
            
            return ChernSimonsNumber;
            
        }
        
        DOUBLE NCS(GaugeLinks *U,GaugeTransformations *G){
            return NCS(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,G);
        }
        
        DOUBLE NCS(GaugeTransformations *G){
            return NCS(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,G);
        }
        
        
        DOUBLE NCS(){
            return NCS(GaugeTransformation::G);
        }
        
    }
    
}
