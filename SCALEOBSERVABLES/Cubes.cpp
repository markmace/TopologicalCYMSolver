namespace ScaleObservables{
    
    // NUMBER OF BASES //
    const INT NumberOfBases=30;
    
    // CUBES SET UP //
    INT SETUP_CUBES=0;
    
    namespace Cubes{
        
        // BASIS DEFINED BY 3 VECTORS //
        struct Basis{
            INT v[3][3];
        };
        
        // SET OF BASES //
        Basis BasisSet[NumberOfBases];
        
        // GENERATE PERMUTATION OF BASIS //
        void GeneratePermutation(INT ux,INT uy,INT uz,INT vx,INT vy,INT vz,INT wx,INT wy,INT wz,Basis &B0,Basis &B1,Basis &B2,Basis &B3,Basis &B4,Basis &B5){
            
            B0.v[0][0]=ux;  B0.v[1][0]=vx;  B0.v[2][0]=wx;
            B0.v[0][1]=uy;  B0.v[1][1]=vy;  B0.v[2][1]=wy;
            B0.v[0][2]=uz;  B0.v[1][2]=vz;  B0.v[2][2]=wz;
            
            B1.v[0][0]=ux;  B1.v[1][0]=vx;  B1.v[2][0]=wx;
            B1.v[0][2]=uy;  B1.v[1][2]=vy;  B1.v[2][2]=wy;
            B1.v[0][1]=uz;  B1.v[1][1]=vz;  B1.v[2][1]=wz;
            
            B2.v[0][1]=ux;  B2.v[1][1]=vx;  B2.v[2][1]=wx;
            B2.v[0][0]=uy;  B2.v[1][0]=vy;  B2.v[2][0]=wy;
            B2.v[0][2]=uz;  B2.v[1][2]=vz;  B2.v[2][2]=wz;
            
            B3.v[0][1]=ux;  B3.v[1][1]=vx;  B3.v[2][1]=wx;
            B3.v[0][2]=uy;  B3.v[1][2]=vy;  B3.v[2][2]=wy;
            B3.v[0][0]=uz;  B3.v[1][0]=vz;  B3.v[2][0]=wz;
            
            B4.v[0][2]=ux;  B4.v[1][2]=vx;  B4.v[2][2]=wx;
            B4.v[0][0]=uy;  B4.v[1][0]=vy;  B4.v[2][0]=wy;
            B4.v[0][1]=uz;  B4.v[1][1]=vz;  B4.v[2][1]=wz;
            
            B5.v[0][2]=ux;  B5.v[1][2]=vx;  B5.v[2][2]=wx;
            B5.v[0][1]=uy;  B5.v[1][1]=vy;  B5.v[2][1]=wy;
            B5.v[0][0]=uz;  B5.v[1][0]=vz;  B5.v[2][0]=wz;
            
        }
        
        
        ////////////////
        // BASIS SETS //
        ////////////////
        
        void Setup(){
            
            // BASIS VECTORS //
            INT ux,uy,uz,vx,vy,vz,wx,wy,wz;
            
            // SET 1 VECTORS //
            ux=3; vx=-4; wx=0;
            uy=4; vy=3;  wy=0;
            uz=0; vz=0;  wz=5;
            
            GeneratePermutation(ux,uy,uz,vx,vy,vz,wx,wy,wz,BasisSet[0],BasisSet[1],BasisSet[2],BasisSet[3],BasisSet[4],BasisSet[5]);
            
            // SET 2 VECTORS //
            ux=5;  vx=-12; wx=0;
            uy=12; vy=5;   wy=0;
            uz=0;  vz=0;   wz=13;
            
            GeneratePermutation(ux,uy,uz,vx,vy,vz,wx,wy,wz,BasisSet[6],BasisSet[7],BasisSet[8],BasisSet[9],BasisSet[10],BasisSet[11]);
            
            // SET 3 VECTORS //
            ux=2; vx=3;  wx=6;
            uy=3; vy=-6; wy=2;
            uz=6; vz=2;  wz=-3;
            
            GeneratePermutation(ux,uy,uz,vx,vy,vz,wx,wy,wz,BasisSet[12],BasisSet[13],BasisSet[14],BasisSet[15],BasisSet[16],BasisSet[17]);
            
            // SET 4 VECTORS //
            ux=1; vx=8;   wx=4;
            uy=4; vy=-4;  wy=7;
            uz=8; vz=1;   wz=-4;
            
            GeneratePermutation(ux,uy,uz,vx,vy,vz,wx,wy,wz,BasisSet[18],BasisSet[19],BasisSet[20],BasisSet[21],BasisSet[22],BasisSet[23]);
            
            // SET 5 VECTORS //
            ux=1; vx=0;  wx=0;
            uy=0; vy=1;  wy=0;
            uz=0; vz=0;  wz=1;
            
            GeneratePermutation(ux,uy,uz,vx,vy,vz,wx,wy,wz,BasisSet[24],BasisSet[25],BasisSet[26],BasisSet[27],BasisSet[28],BasisSet[29]);
            
            // SETUP FLAG //
            SETUP_CUBES=1;
            
        }
        
    }
    
    void GetSquare(INT &x1,INT &y1,INT &z1,INT &x2,INT &y2,INT &z2,INT &x3,INT &y3,INT &z3,INT MAX_SIDE_LENGTH_SQUARED){
        
        // BASIS AND VECTOR SET //
        INT bI,uI,vI;
        
        // BASIS VECTOR NORM AND DOT PRODUCT //
        INT uLengthSqr,vLengthSqr,uDotv;
        
        // SELECT A RANDOM BASE POINT //
        x1=INT(RandomNumberGenerator::rng()*Lattice::N[0]);
        y1=INT(RandomNumberGenerator::rng()*Lattice::N[1]);
        z1=INT(RandomNumberGenerator::rng()*Lattice::N[2]);
        
        
        // SETUP //
        if(SETUP_CUBES==0){
            Cubes::Setup();
        }
        
        
        // DETERMINE A BASIS TO USE  //
        INT CHECK_LENGTH_FLAG=0;
        
        while(CHECK_LENGTH_FLAG==0){
            

            // CHOOSE A RANDOM BASIS //
            bI=INT(NumberOfBases*RandomNumberGenerator::rng());
            
            // CHOOSE BASE VECTORS //
            uI=0; vI=0;
            
            while(uI==vI){
                uI=INT(3*RandomNumberGenerator::rng()); vI=INT(3*RandomNumberGenerator::rng());
            }
            
            // GET LENGTH OF BASIS VECTORS //
            uLengthSqr=SQR(Cubes::BasisSet[bI].v[uI][0])+SQR(Cubes::BasisSet[bI].v[uI][1])+SQR(Cubes::BasisSet[bI].v[uI][2]);
            vLengthSqr=SQR(Cubes::BasisSet[bI].v[vI][0])+SQR(Cubes::BasisSet[bI].v[vI][1])+SQR(Cubes::BasisSet[bI].v[vI][2]);
            
            // COMPUTE SCALAR PRODUCT //
            uDotv=Cubes::BasisSet[bI].v[uI][0]*Cubes::BasisSet[bI].v[vI][0]+Cubes::BasisSet[bI].v[uI][1]*Cubes::BasisSet[bI].v[vI][1]+Cubes::BasisSet[bI].v[uI][2]*Cubes::BasisSet[bI].v[vI][2];
            
            // CHECK ORTHOGONALITY AND NORMALIZATION //
            if(uDotv!=0 || uLengthSqr!=vLengthSqr){
                std::cerr << "#VECTORS INCORRECTLY ASSIGNED -- FAILED ORTHOGONALITY AND/OR EQUAL MAGNITUDE TEST" << std::endl;
                exit(0);
            }

            // CHECK THAT LENGTH IS ACCEPTABLE -- IF NOT REPEAT //
            if(uLengthSqr<MAX_SIDE_LENGTH_SQUARED && vLengthSqr<MAX_SIDE_LENGTH_SQUARED){
                CHECK_LENGTH_FLAG=1;
            }
        
        }
    
        
        // DETERMINE SCALING FACTORS //
        INT MaxScalingFactor=INT(sqrt(MAX_SIDE_LENGTH_SQUARED/uLengthSqr));
        
        // SCALING FACTORS //
        INT alpha,beta;
        
        // SAMPLE SQUARE WITH ACCEPTABLE SIDE LENGTH //
        INT CHECK_ACCEPTANCE=0;
        
        while(CHECK_ACCEPTANCE==0){
            
            alpha=INT((2*MaxScalingFactor+1)*RandomNumberGenerator::rng())-MaxScalingFactor;
            beta =INT((2*MaxScalingFactor+1)*RandomNumberGenerator::rng())-MaxScalingFactor;
            
            // CHECK LENGTH //
            if(((SQR(alpha)*uLengthSqr+SQR(beta)*vLengthSqr)<=MAX_SIDE_LENGTH_SQUARED) && ((SQR(alpha)*uLengthSqr+SQR(beta)*vLengthSqr)>0)){
                CHECK_ACCEPTANCE=1;
            }
            
        }
        
        ///////////////////////////////////
        // GENERATE DISPLACEMENT VECTORS //
        ///////////////////////////////////
        
        // SET SECOND POINT //
        
        x2=MOD(x1+alpha*Cubes::BasisSet[bI].v[uI][0]+beta*Cubes::BasisSet[bI].v[vI][0],Lattice::N[0]);
        y2=MOD(y1+alpha*Cubes::BasisSet[bI].v[uI][1]+beta*Cubes::BasisSet[bI].v[vI][1],Lattice::N[1]);
        z2=MOD(z1+alpha*Cubes::BasisSet[bI].v[uI][2]+beta*Cubes::BasisSet[bI].v[vI][2],Lattice::N[2]);
        
        // SET THIRD POINT //
        
        x3=MOD(x2-beta*Cubes::BasisSet[bI].v[uI][0]+alpha*Cubes::BasisSet[bI].v[vI][0],Lattice::N[0]);
        y3=MOD(y2-beta*Cubes::BasisSet[bI].v[uI][1]+alpha*Cubes::BasisSet[bI].v[vI][1],Lattice::N[1]);
        z3=MOD(z2-beta*Cubes::BasisSet[bI].v[uI][2]+alpha*Cubes::BasisSet[bI].v[vI][2],Lattice::N[2]);
                
    }
}