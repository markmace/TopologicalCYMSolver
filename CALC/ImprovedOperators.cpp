namespace ImprovedOperators{
    
    // COMPUTE PARALLEL TRANSPORT OF SU(N) MATRIX //
    void ParallelTransport(SU_Nc_FUNDAMENTAL_FORMAT *U,SU_Nc_FUNDAMENTAL_FORMAT *Q){
        
        SUNcGroup::AdvancedOperations::UUD(U,Q,U,Q);
    }
    
    DOUBLE Compute(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        // tr E.B //
        DOUBLE EDotB=0.0;

        #pragma omp parallel
        {
            
            // AVERAGE PLAQUETTES //
            SU_Nc_FUNDAMENTAL_FORMAT UxyAvg[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UyzAvg[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UzxAvg[SUNcGroup::MatrixSize];
            
            // ELECTRIC AND MAGNETIC FIELDS //
            DOUBLE Ex[SUNcAlgebra::VectorSize]; DOUBLE Ey[SUNcAlgebra::VectorSize]; DOUBLE Ez[SUNcAlgebra::VectorSize];
            DOUBLE Bx[SUNcAlgebra::VectorSize]; DOUBLE By[SUNcAlgebra::VectorSize]; DOUBLE Bz[SUNcAlgebra::VectorSize];
            
            //////////////////////
            // ALLOCATE BUFFERS //
            //////////////////////
            
            // ELEMENTARY PLAQUETTES //
            
            // Ux+1,y+1 Uy+1,z+1 Uz+1,x+1
            SU_Nc_FUNDAMENTAL_FORMAT UBoxXY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxYZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxZX[SUNcGroup::MatrixSize];
            
            // Ux-1,y+1 Uy-1,z+1 Uz-1,x+1
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMXY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMYZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMZX[SUNcGroup::MatrixSize];
            
            // Ux-1,y-1 Uy-1,z-1 Uz-1,x-1
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMXMY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMYMZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxMZMX[SUNcGroup::MatrixSize];
            
            // Ux+1,y-1 Uy+1,z-1 Uz+1,x-1
            SU_Nc_FUNDAMENTAL_FORMAT UBoxXMY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxYMZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UBoxZMX[SUNcGroup::MatrixSize];
            
            // RECTANGULAR PLAQUETTES //
            
            // Ux+2,y+1 Uy+2,z+1 Uz+2,y+1
            SU_Nc_FUNDAMENTAL_FORMAT URectX2Y[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectY2Z[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectZ2X[SUNcGroup::MatrixSize];
            
            // Ux+1,y+2 Uy+1,z+2 Uz+1,x+2
            SU_Nc_FUNDAMENTAL_FORMAT URectXY2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectYZ2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectZX2[SUNcGroup::MatrixSize];
            
            // Ux-2,y+1 Uy-2,z+1 Uz-2,x+1
            SU_Nc_FUNDAMENTAL_FORMAT URectMX2Y[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMY2Z[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMZ2X[SUNcGroup::MatrixSize];
            
            // Ux-1,y+2 Uy-1,z+2 Uz-1,x+2
            SU_Nc_FUNDAMENTAL_FORMAT URectMXY2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMYZ2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMZX2[SUNcGroup::MatrixSize];
            
            // Ux-2,y-1 Uy-2,z-1 Uz-2,x-1
            SU_Nc_FUNDAMENTAL_FORMAT URectMX2MY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMY2MZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMZ2MX[SUNcGroup::MatrixSize];
            
            // Ux-1,y-2 Uy-1,z-2 Uz-1,x-2
            SU_Nc_FUNDAMENTAL_FORMAT URectMXMY2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMYMZ2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectMZMX2[SUNcGroup::MatrixSize];
            
            // Ux+2,y-1 Uy+2,z-1 Uz+2,y-1
            SU_Nc_FUNDAMENTAL_FORMAT URectX2MY[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectY2MZ[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectZ2MX[SUNcGroup::MatrixSize];
            
            // Ux+1,y-2 Uy+1,z-2 Uz+1,x-2
            SU_Nc_FUNDAMENTAL_FORMAT URectXMY2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectYMZ2[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT URectZMX2[SUNcGroup::MatrixSize];
            
            
            // PARALLEL TRANSPORTERS //
            
            // U(x+1)->(x) //
            SU_Nc_ADJOINT_FORMAT UAdjX[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjY[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjZ[SUNcGroup::AdjointSize];
            
            // U(x-1)->(x) //
            SU_Nc_FUNDAMENTAL_FORMAT UxM[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UyM[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT UzM[SUNcGroup::MatrixSize];
            
            SU_Nc_ADJOINT_FORMAT UAdjxM[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjyM[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjzM[SUNcGroup::AdjointSize];
            
            // U(x-2)->(x) //
            SU_Nc_FUNDAMENTAL_FORMAT Ux2M[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT Uy2M[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT Uz2M[SUNcGroup::MatrixSize];
            
            SU_Nc_ADJOINT_FORMAT UAdjx2M[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjy2M[SUNcGroup::AdjointSize];
            SU_Nc_ADJOINT_FORMAT UAdjz2M[SUNcGroup::AdjointSize];
            
            // COMPUTE O(a^2) IMPROVED MAGNETIC AND ELECTRIC FIELDS //
            #pragma omp for reduction( + : EDotB)
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        /////////////////////////////////
                        //  GET PARALLEL TRANSPORTERS  //
                        /////////////////////////////////
                        
                        // SINGLE TRANSPORT (x+1)->x //
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,0),UAdjX);
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,1),UAdjY);
                        SUNcGroup::Operations::GetAdjoint(U->Get(x,y,z,2),UAdjZ);
                        
                        
                        // SINGLE TRANSPORT (x-1)->x //
                        SUNcGroup::Operations::Inverse(U->Get(x-1,y,z,0),UxM);
                        SUNcGroup::Operations::Inverse(U->Get(x,y-1,z,1),UyM);
                        SUNcGroup::Operations::Inverse(U->Get(x,y,z-1,2),UzM);
                        
                        SUNcGroup::Operations::GetAdjoint(UxM,UAdjxM);
                        SUNcGroup::Operations::GetAdjoint(UyM,UAdjyM);
                        SUNcGroup::Operations::GetAdjoint(UzM,UAdjzM);
                        
                        // DOUBLE TRANSPORT (x-2) -> x //
                        SUNcGroup::Operations::DD(U->Get(x-1,y,z,0),U->Get(x-2,y,z,0),Ux2M);
                        SUNcGroup::Operations::DD(U->Get(x,y-1,z,1),U->Get(x,y-2,z,1),Uy2M);
                        SUNcGroup::Operations::DD(U->Get(x,y,z-1,2),U->Get(x,y,z-2,2),Uz2M);
                        
                        
                        SUNcGroup::Operations::GetAdjoint(Ux2M,UAdjx2M);
                        SUNcGroup::Operations::GetAdjoint(Uy2M,UAdjy2M);
                        SUNcGroup::Operations::GetAdjoint(Uz2M,UAdjz2M);
                        
                        
                        ///////////////////////////////////
                        // COMPUTE ELEMENTARY PLAQUETTES //
                        ///////////////////////////////////
                        
                        // UPPER RIGHT BOX //
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,0),U->Get(x+1,y,z,1),U->Get(x,y+1,z,0),U->Get(x,y,z,1),UBoxXY);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,1),U->Get(x,y+1,z,2),U->Get(x,y,z+1,1),U->Get(x,y,z,2),UBoxYZ);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z,2),U->Get(x,y,z+1,0),U->Get(x+1,y,z,2),U->Get(x,y,z,0),UBoxZX);
                        
                        // UPPER LEFT BOX //
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x-1,y,z,0),U->Get(x,y,z,1),U->Get(x-1,y+1,z,0),U->Get(x-1,y,z,1),UBoxMXY);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y-1,z,1),U->Get(x,y,z,2),U->Get(x,y-1,z+1,1),U->Get(x,y-1,z,2),UBoxMYZ);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z-1,2),U->Get(x,y,z,0),U->Get(x+1,y,z-1,2),U->Get(x,y,z-1,0),UBoxMZX);
                        
                        ParallelTransport(UxM,UBoxMXY);
                        ParallelTransport(UyM,UBoxMYZ);
                        ParallelTransport(UzM,UBoxMZX);
                        
                        // LOWER LEFT BOX //
                        SUNcGroup::AdvancedOperations::DDUU(U->Get(x-1,y,z,0),U->Get(x-1,y-1,z,1),U->Get(x-1,y-1,z,0),U->Get(x,y-1,z,1),UBoxMXMY);
                        SUNcGroup::AdvancedOperations::DDUU(U->Get(x,y-1,z,1),U->Get(x,y-1,z-1,2),U->Get(x,y-1,z-1,1),U->Get(x,y,z-1,2),UBoxMYMZ);
                        SUNcGroup::AdvancedOperations::DDUU(U->Get(x,y,z-1,2),U->Get(x-1,y,z-1,0),U->Get(x-1,y,z-1,2),U->Get(x-1,y,z,0),UBoxMZMX);
                        
                        // LOWER RIGHT BOX //
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y-1,z,0),U->Get(x+1,y-1,z,1),U->Get(x,y,z,0),U->Get(x,y-1,z,1),UBoxXMY);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x,y,z-1,1),U->Get(x,y+1,z-1,2),U->Get(x,y,z,1),U->Get(x,y,z-1,2),UBoxYMZ);
                        SUNcGroup::AdvancedOperations::UUDD(U->Get(x-1,y,z,2),U->Get(x-1,y,z+1,0),U->Get(x,y,z,2),U->Get(x-1,y,z,0),UBoxZMX);
                        
                        ParallelTransport(UyM,UBoxXMY);
                        ParallelTransport(UzM,UBoxYMZ);
                        ParallelTransport(UxM,UBoxZMX);
                        
                        ////////////////////////////////////
                        // COMPUTE RECTANGULAR PLAQUETTES //
                        ////////////////////////////////////
                        
                        // UPPER RIGHT RECTANGLES //
                        // 1 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,0),U->Get(x+1,y,z,0),U->Get(x+2,y,z,1),U->Get(x+1,y+1,z,0),U->Get(x,y+1,z,0),U->Get(x,y,z,1),URectX2Y);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,1),U->Get(x,y+1,z,1),U->Get(x,y+2,z,2),U->Get(x,y+1,z+1,1),U->Get(x,y,z+1,1),U->Get(x,y,z,2),URectY2Z);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,2),U->Get(x,y,z+1,2),U->Get(x,y,z+2,0),U->Get(x+1,y,z+1,2),U->Get(x+1,y,z,2),U->Get(x,y,z,0),URectZ2X);
                        // 8 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,0),U->Get(x+1,y,z,1),U->Get(x+1,y+1,z,1),U->Get(x,y+2,z,0),U->Get(x,y+1,z,1),U->Get(x,y,z,1),URectXY2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,1),U->Get(x,y+1,z,2),U->Get(x,y+1,z+1,2),U->Get(x,y,z+2,1),U->Get(x,y,z+1,2),U->Get(x,y,z,2),URectYZ2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z,2),U->Get(x,y,z+1,0),U->Get(x+1,y,z+1,0),U->Get(x+2,y,z,2),U->Get(x+1,y,z,0),U->Get(x,y,z,0),URectZX2);
                        
                        // UPPER LEFT RECTANGLES //
                        // 7 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x-1,y,z,0),U->Get(x,y,z,1),U->Get(x,y+1,z,1),U->Get(x-1,y+2,z,0),U->Get(x-1,y+1,z,1),U->Get(x-1,y,z,1),URectMXY2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y-1,z,1),U->Get(x,y,z,2),U->Get(x,y,z+1,2),U->Get(x,y-1,z+2,1),U->Get(x,y-1,z+1,2),U->Get(x,y-1,z,2),URectMYZ2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z-1,2),U->Get(x,y,z,0),U->Get(x+1,y,z,0),U->Get(x+2,y,z-1,2),U->Get(x+1,y,z-1,0),U->Get(x,y,z-1,0),URectMZX2);
                        
                        ParallelTransport(UxM,URectMXY2);
                        ParallelTransport(UyM,URectMYZ2);
                        ParallelTransport(UzM,URectMZX2);
                        
                        // 6 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x-2,y,z,0),U->Get(x-1,y,z,0),U->Get(x,y,z,1),U->Get(x-1,y+1,z,0),U->Get(x-2,y+1,z,0),U->Get(x-2,y,z,1),URectMX2Y);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y-2,z,1),U->Get(x,y-1,z,1),U->Get(x,y,z,2),U->Get(x,y-1,z+1,1),U->Get(x,y-2,z+1,1),U->Get(x,y-2,z,2),URectMY2Z);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z-2,2),U->Get(x,y,z-1,2),U->Get(x,y,z,0),U->Get(x+1,y,z-1,2),U->Get(x+1,y,z-2,2),U->Get(x,y,z-2,0),URectMZ2X);
                        
                        ParallelTransport(Ux2M,URectMX2Y);
                        ParallelTransport(Uy2M,URectMY2Z);
                        ParallelTransport(Uz2M,URectMZ2X);
                        
                        
                        // LOWER LEFT RECTANGLES //
                        // 4 //
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x-1,y,z,0),U->Get(x-1,y-1,z,1),U->Get(x-1,y-2,z,1),U->Get(x-1,y-2,z,0),U->Get(x,y-2,z,1),U->Get(x,y-1,z,1),URectMXMY2);
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x,y-1,z,1),U->Get(x,y-1,z-1,2),U->Get(x,y-1,z-2,2),U->Get(x,y-1,z-2,1),U->Get(x,y,z-2,2),U->Get(x,y,z-1,2),URectMYMZ2);
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x,y,z-1,2),U->Get(x-1,y,z-1,0),U->Get(x-2,y,z-1,0),U->Get(x-2,y,z-1,2),U->Get(x-2,y,z,0),U->Get(x-1,y,z,0),URectMZMX2);
                        // 5 //
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x-1,y,z,0),U->Get(x-2,y,z,0),U->Get(x-2,y-1,z,1),U->Get(x-2,y-1,z,0),U->Get(x-1,y-1,z,0),U->Get(x,y-1,z,1),URectMX2MY);
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x,y-1,z,1),U->Get(x,y-2,z,1),U->Get(x,y-2,z-1,2),U->Get(x,y-2,z-1,1),U->Get(x,y-1,z-1,1),U->Get(x,y,z-1,2),URectMY2MZ);
                        SUNcGroup::AdvancedOperations::DDDUUU(U->Get(x,y,z-1,2),U->Get(x,y,z-2,2),U->Get(x-1,y,z-2,0),U->Get(x-1,y,z-2,2),U->Get(x-1,y,z-1,2),U->Get(x-1,y,z,0),URectMZ2MX);
                        
                        
                        // LOWER RIGHT RECTANGLES //
                        // 2 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y-1,z,0),U->Get(x+1,y-1,z,0),U->Get(x+2,y-1,z,1),U->Get(x+1,y,z,0),U->Get(x,y,z,0),U->Get(x,y-1,z,1),URectX2MY);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z-1,1),U->Get(x,y+1,z-1,1),U->Get(x,y+2,z-1,2),U->Get(x,y+1,z,1),U->Get(x,y,z,1),U->Get(x,y,z-1,2),URectY2MZ);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x-1,y,z,2),U->Get(x-1,y,z+1,2),U->Get(x-1,y,z+2,0),U->Get(x,y,z+1,2),U->Get(x,y,z,2),U->Get(x-1,y,z,0),URectZ2MX);
                        
                        ParallelTransport(UyM,URectX2MY);
                        ParallelTransport(UzM,URectY2MZ);
                        ParallelTransport(UxM,URectZ2MX);
                        // 3 //
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y-2,z,0),U->Get(x+1,y-2,z,1),U->Get(x+1,y-1,z,1),U->Get(x,y,z,0),U->Get(x,y-1,z,1),U->Get(x,y-2,z,1),URectXMY2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x,y,z-2,1),U->Get(x,y+1,z-2,2),U->Get(x,y+1,z-1,2),U->Get(x,y,z,1),U->Get(x,y,z-1,2),U->Get(x,y,z-2,2),URectYMZ2);
                        SUNcGroup::AdvancedOperations::UUUDDD(U->Get(x-2,y,z,2),U->Get(x-2,y,z+1,0),U->Get(x-1,y,z+1,0),U->Get(x,y,z,2),U->Get(x-1,y,z,0),U->Get(x-2,y,z,0),URectZMX2);
                        
                        ParallelTransport(Uy2M,URectXMY2);
                        ParallelTransport(Uz2M,URectYMZ2);
                        ParallelTransport(Ux2M,URectZMX2);
                        
                        
                        //////////////////////////////////////////////////
                        // COMPUTE AVERAGE ELECTRIC AND MAGNETIC FIELDS //
                        //////////////////////////////////////////////////
                        
                        // COMPUTE ELECTRIC FIELDS //
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            Ex[a]=0.0; Ey[a]=0.0; Ez[a]=0.0;
                            
                            for(INT b=0;b<SUNcAlgebra::VectorSize;b++){
                                
                                Ex[a]+=-(1.0/12.0)*UAdjX[SUNcGroup::AdjIndex(a,b)]*E->Get(x+1,y,z,0,b)[0] + (7.0/12.0)*DELTA(a,b)*E->Get(x,y,z,0,a)[0] + (7.0/12.0)*UAdjxM[SUNcGroup::AdjIndex(a,b)]*E->Get(x-1,y,z,0,b)[0] - (1.0/12.0)*UAdjx2M[SUNcGroup::AdjIndex(a,b)]*E->Get(x-2,y,z,0,b)[0];
                                
                                Ey[a]+=-(1.0/12.0)*UAdjY[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y+1,z,1,b)[0] + (7.0/12.0)*DELTA(a,b)*E->Get(x,y,z,1,a)[0] + (7.0/12.0)*UAdjyM[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y-1,z,1,b)[0] - (1.0/12.0)*UAdjy2M[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y-2,z,1,b)[0];
                                
                                Ez[a]+=-(1.0/12.0)*UAdjZ[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y,z+1,2,b)[0] + (7.0/12.0)*DELTA(a,b)*E->Get(x,y,z,2,a)[0] + (7.0/12.0)*UAdjzM[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y,z-1,2,b)[0] - (1.0/12.0)*UAdjz2M[SUNcGroup::AdjIndex(a,b)]*E->Get(x,y,z-2,2,b)[0];
                                
                            }
                            
                        }
                        
                        // COMPUTE MAGNETIC FIELDS //
                        for(INT s=0;s<SUNcGroup::MatrixSize;s++){
                            
                            UxyAvg[s]=(5.0/12.0)*(UBoxXY[s]+UBoxMXY[s]+UBoxMXMY[s]+UBoxXMY[s])-(1.0/24.0)*(URectX2Y[s]+URectXY2[s]+URectMX2Y[s]+URectMXY2[s]+URectMX2MY[s]+URectMXMY2[s]+URectXMY2[s]+URectX2MY[s]);
                            UyzAvg[s]=(5.0/12.0)*(UBoxYZ[s]+UBoxMYZ[s]+UBoxMYMZ[s]+UBoxYMZ[s])-(1.0/24.0)*(URectY2Z[s]+URectYZ2[s]+URectMY2Z[s]+URectMYZ2[s]+URectMY2MZ[s]+URectMYMZ2[s]+URectYMZ2[s]+URectY2MZ[s]);
                            UzxAvg[s]=(5.0/12.0)*(UBoxZX[s]+UBoxMZX[s]+UBoxMZMX[s]+UBoxZMX[s])-(1.0/24.0)*(URectZ2X[s]+URectZX2[s]+URectMZ2X[s]+URectMZX2[s]+URectMZ2MX[s]+URectMZMX2[s]+URectZMX2[s]+URectZ2MX[s]);
                            
                        }
                        
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UyzAvg,Bx);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UzxAvg,By);
                        SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0),UxyAvg,Bz);
                        
                        
                        // COMPUTE TRACE OF E.B //
                        
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            EDotB+=Ex[a]*Bx[a]+Ey[a]*By[a]+Ez[a]*Bz[a];
                        }
                        
                    }
                    
                }
            }
            
        } // END PARALLEL //
        
        return EDotB;
        
        
    }
    
    DOUBLE Compute(GaugeLinks *U,ElectricFields *E){
        
        return Compute(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,U,E);
        
    }
    
    DOUBLE Compute(){
        
        return Compute(GLinks::U,EFields::E);
        
    }

}