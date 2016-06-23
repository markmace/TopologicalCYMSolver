namespace ScaleObservables{
    
    //////////////////////////////////
    // SINGLE GAUGE LINK OPERATIONS //
    //////////////////////////////////
    
    void GoXPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UU(Ux,U->Get(x+Var,y,z,0),UxTmp);
            }
            else{
                SUNcGroup::Operations::UU(UxTmp,U->Get(x+Var,y,z,0),Ux);
            }
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
            
        }
        
    }
    
    void GoXMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UD(Ux,U->Get(x-Var-1,y,z,0),UxTmp);
            }
            else{
                SUNcGroup::Operations::UD(UxTmp,U->Get(x-Var-1,y,z,0),Ux);
            }
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
        }
        
    }
    
    void GoYPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UU(Ux,U->Get(x,y+Var,z,1),UxTmp);
            }
            else{
                SUNcGroup::Operations::UU(UxTmp,U->Get(x,y+Var,z,1),Ux);
            }
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
        }
        
    }
    
    void GoYMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UD(Ux,U->Get(x,y-Var-1,z,1),UxTmp);
            }
            else{
                SUNcGroup::Operations::UD(UxTmp,U->Get(x,y-Var-1,z,1),Ux);
            }
            
            
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
        }
        
    }
    
    void GoZPlus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UU(Ux,U->Get(x,y,z+Var,2),UxTmp);
            }
            else{
                SUNcGroup::Operations::UU(UxTmp,U->Get(x,y,z+Var,2),Ux);
            }
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
        }
        
    }
    
    void GoZMinus(INT x,INT y,INT z,INT L,SU_Nc_FUNDAMENTAL_FORMAT *Ux,GaugeLinks *U){
        
        SU_Nc_FUNDAMENTAL_FORMAT UxTmp[SUNcGroup::MatrixSize];
        
        for(INT Var=0;Var<L;Var++){
            
            if(Var%2==0){
                SUNcGroup::Operations::UD(Ux,U->Get(x,y,z-Var-1,2),UxTmp);
            }
            else{
                SUNcGroup::Operations::UD(UxTmp,U->Get(x,y,z-Var-1,2),Ux);
            }
        }
        
        if((L+1)%2==0){
            COPY_SUNcMatrix(Ux,UxTmp);
        }
        
    }

    //////////////////////////////////
    // DISTANCE MEASURES AND PATHES //
    //////////////////////////////////
    
    // GET PERIODIC DISTANCE IN ONE DIMENSION //
    
    DOUBLE GetDistance(INT x1,INT x2,INT N){
        
        INT dX=x2-x1;
        
        INT dXLat=( dX + N*((dX<0)?(1):(0)));
        
        return ((dXLat<=(N/2)) ? (dXLat) : -(N-dXLat));
        
    }

    
    // COMPUTE SHORTEST SQUARE DISTANCE BETWEEN TWO POINTS ON PERIODIC LATTICE //
    DOUBLE GetSqrDistance(INT x1,INT y1,INT z1,INT x2,INT y2,INT z2){
        
        DOUBLE dX=GetDistance(x1,x2,Lattice::N[0]);
        DOUBLE dY=GetDistance(y1,y2,Lattice::N[1]);
        DOUBLE dZ=GetDistance(z1,z2,Lattice::N[2]);
        
        return SQR(dX)*SQR(Lattice::a[0])+SQR(dY)*SQR(Lattice::a[1])+SQR(dZ)*SQR(Lattice::a[2]);
        
    }
    
    // COMPUTE SHORTEST DISTANCE BETWEEN TWO POINTS ON PERIODIC LATTICE //
    DOUBLE GetDistance(INT x1,INT y1,INT z1,INT x2,INT y2,INT z2){
        
        return sqrt(GetSqrDistance(x1,y1,z1,x2,y2,z2));
        
    }
    
    // GET STEP ALONG PATH CLOSEST TO GEODESIC //
    // (-1,0,0) -> -1 | (+1,0,0) -> +1 | (0,-1,0) -> -2 | (0,+1,0) ->+2 | (0,0,-1) -> -3 | (0,0,+1) -> +3 //
    
    INT GetStep(INT xStart,INT yStart,INT zStart,INT xEnd,INT yEnd,INT zEnd,INT xPos,INT yPos,INT zPos){
        
        // MINIMAL PATH LENGTH //
        DOUBLE Des=GetSqrDistance(xStart,yStart,zStart,xEnd,yEnd,zEnd);
        
        // ORIGINAL DISTANCES TO START AND END POINT //
        DOUBLE DsOriginal=GetSqrDistance(xStart,yStart,zStart,xPos,yPos,zPos);
        DOUBLE DeOriginal=GetSqrDistance(xEnd,yEnd,zEnd,xPos,yPos,zPos);
        
        //////////////////////////////
        // CHECK ALL POSSIBLE STEPS //
        //////////////////////////////
        
        // CHECK DISTANCE TO START POINT //
        DOUBLE DsXp=GetSqrDistance(xStart,yStart,zStart,MOD((xPos+1),Lattice::N[0]),yPos,zPos);
        DOUBLE DsXm=GetSqrDistance(xStart,yStart,zStart,MOD((xPos-1),Lattice::N[0]),yPos,zPos);
        
        DOUBLE DsYp=GetSqrDistance(xStart,yStart,zStart,xPos,MOD((yPos+1),Lattice::N[1]),zPos);
        DOUBLE DsYm=GetSqrDistance(xStart,yStart,zStart,xPos,MOD((yPos-1),Lattice::N[1]),zPos);
        
        DOUBLE DsZp=GetSqrDistance(xStart,yStart,zStart,xPos,yPos,MOD((zPos+1),Lattice::N[2]));
        DOUBLE DsZm=GetSqrDistance(xStart,yStart,zStart,xPos,yPos,MOD((zPos-1),Lattice::N[2]));
        
        // CHECK DISTANCE TO END POINT //
        DOUBLE DeXp=GetSqrDistance(xEnd,yEnd,zEnd,MOD((xPos+1),Lattice::N[0]),yPos,zPos);
        DOUBLE DeXm=GetSqrDistance(xEnd,yEnd,zEnd,MOD((xPos-1),Lattice::N[0]),yPos,zPos);
        
        DOUBLE DeYp=GetSqrDistance(xEnd,yEnd,zEnd,xPos,MOD((yPos+1),Lattice::N[1]),zPos);
        DOUBLE DeYm=GetSqrDistance(xEnd,yEnd,zEnd,xPos,MOD((yPos-1),Lattice::N[1]),zPos);
        
        DOUBLE DeZp=GetSqrDistance(xEnd,yEnd,zEnd,xPos,yPos,MOD((zPos+1),Lattice::N[2]));
        DOUBLE DeZm=GetSqrDistance(xEnd,yEnd,zEnd,xPos,yPos,MOD((zPos-1),Lattice::N[2]));
        
        
        // GET DISTANCE FROM IDEAL PATH //
        DOUBLE hXp=0.5*(DsXp+DeXp)-0.25*SQR(DsXp-DeXp)/Des-0.25*Des;
        DOUBLE hXm=0.5*(DsXm+DeXm)-0.25*SQR(DsXm-DeXm)/Des-0.25*Des;
        
        DOUBLE hYp=0.5*(DsYp+DeYp)-0.25*SQR(DsYp-DeYp)/Des-0.25*Des;
        DOUBLE hYm=0.5*(DsYm+DeYm)-0.25*SQR(DsYm-DeYm)/Des-0.25*Des;
        
        DOUBLE hZp=0.5*(DsZp+DeZp)-0.25*SQR(DsZp-DeZp)/Des-0.25*Des;
        DOUBLE hZm=0.5*(DsZm+DeZm)-0.25*SQR(DsZm-DeZm)/Des-0.25*Des;
        
        
        // DETERMINE MINIMUM DISTANCE FROM ORIGINAL PATH //
        DOUBLE hMin=SQR(Lattice::N[0]*Lattice::a[0])+SQR(Lattice::N[1]*Lattice::a[1])+SQR(Lattice::N[2]*Lattice::a[2]);
        
        if(DeXp<=DeOriginal){
            hMin=std::min(hMin,hXp);
        }
        if(DeXm<=DeOriginal){
            hMin=std::min(hMin,hXm);
        }
        if(DeYp<=DeOriginal){
            hMin=std::min(hMin,hYp);
        }
        if(DeYm<=DeOriginal){
            hMin=std::min(hMin,hYm);
        }
        if(DeZp<=DeOriginal){
            hMin=std::min(hMin,hZp);
        }
        if(DeZm<=DeOriginal){
            hMin=std::min(hMin,hZm);
        }
        
        // CHECK WHICH OPTIONS ARE AVAILABLE //
        unsigned int NumberOfOptions=0;
        
        std::vector<INT> Options;
        
        if(DeXp<=DeOriginal && hXp==hMin){
            NumberOfOptions++; Options.push_back(+1);
        }
        if(DeXm<=DeOriginal && hXm==hMin){
            NumberOfOptions++; Options.push_back(-1);
        }
        if(DeYp<=DeOriginal && hYp==hMin){
            NumberOfOptions++; Options.push_back(+2);
        }
        if(DeYm<=DeOriginal && hYm==hMin){
            NumberOfOptions++; Options.push_back(-2);
        }
        if(DeZp<=DeOriginal && hZp==hMin){
            NumberOfOptions++; Options.push_back(+3);
        }
        if(DeZm<=DeOriginal && hZm==hMin){
            NumberOfOptions++; Options.push_back(-3);
        }
        
        // CHECK THAT OPTIONS ARE AVAILABLE //
        if(NumberOfOptions==0 || NumberOfOptions!=Options.size()){
            
            std::cerr << "#VACATION TIME IS OVER -- NOpt=" << NumberOfOptions << " VS=" << Options.size()  << std::endl;
            
            std::cerr << "#DIST START " << DsXm << " " << DsXp << " " << DsYm << " " << DsYp << " " << DsZm << " " << DsZp << "  ---  " << DsOriginal << std::endl;
            std::cerr << "#DIST END " << DeXm << " " << DeXp << " " << DeYm << " " << DeYp << " " << DeZm << " " << DeZp << "  ---  " << DeOriginal << std::endl;
            std::cerr << "#HEIGHT " << hXm << " " << hXp << " " << hYm << " " << hYp << " " << hZm << " " << hZp << "  ---  " << hMin << std::endl;
            
            exit(0);
        }
        
        // SELECT OPTION //
        INT Choice=INT(NumberOfOptions*RandomNumberGenerator::rng());
        
        return Options[Choice];
        
    }
    
    void GetPath(INT xStart,INT yStart,INT zStart,INT xEnd,INT yEnd,INT zEnd,std::vector<INT>& Path){
        
        // START //
        INT xPos=xStart; INT yPos=yStart; INT zPos=zStart;
        
        // GET INITIAL DISTANCE //
        DOUBLE Distance=GetSqrDistance(xPos,yPos,zPos,xEnd,yEnd,zEnd);
        
        // DETERMINE PATH //
        while(Distance>0){
            
            // GET STEP //
            INT Increment=GetStep(xStart,yStart,zStart,xEnd,yEnd,zEnd,xPos,yPos,zPos);
            
            // MAKE STEP //
            if(Increment==-1){
                xPos=MOD((xPos-1),Lattice::N[0]);
            }
            if(Increment==+1){
                xPos=MOD((xPos+1),Lattice::N[0]);
            }
            if(Increment==-2){
                yPos=MOD((yPos-1),Lattice::N[1]);
            }
            if(Increment==+2){
                yPos=MOD((yPos+1),Lattice::N[1]);
            }
            if(Increment==-3){
                zPos=MOD((zPos-1),Lattice::N[2]);
            }
            if(Increment==+3){
                zPos=MOD((zPos+1),Lattice::N[2]);
            }
            
            // SAVE STEP //
            Path.push_back(Increment);
            
            Distance=GetSqrDistance(xPos,yPos,zPos,xEnd,yEnd,zEnd);
        }
        
        
    }
    
    void GetInversePath(std::vector<INT> OriginalPath, std::vector<INT> &InversePath){
        
        // GET PATH LENGTH //
        INT PathLength=OriginalPath.size();
        
        // MOVE BACKWARDS ALONG ORIGINAL PATH //
        for(INT n=PathLength-1;n>=0;n--){
            InversePath.push_back(-OriginalPath[n]);
        }
        
        // CHECK //
        if(OriginalPath.size()!=InversePath.size()){
            std::cerr << "#ERROR -- INVERSE PATH NOT DETERMINED CORRECTLY" << std::endl;
            exit(0);
        }
        
    }
    
    /////////////////////////////
    // GAUGE LINK TRANSPORTERS //
    /////////////////////////////
    
    
    // DETERMINE TRANSPORTER FOR DEFINED START AND END POINTS //
    void GetTransporter(INT xStart,INT yStart,INT zStart,INT xEnd,INT yEnd,INT zEnd,SU_Nc_FUNDAMENTAL_FORMAT *W,GaugeLinks *U){
        
        // START //
        INT xPos=xStart; INT yPos=yStart; INT zPos=zStart;
        
        // RESET LINK //
        COPY_SUNcMatrix(W,SUNcGroup::UnitMatrix);
        
        // GET INITIAL DISTANCE //
        DOUBLE Distance=GetSqrDistance(xPos,yPos,zPos,xEnd,yEnd,zEnd);
        
        // OUTPUT STARTPOINT OF PATH //
        //std::cout << "#INIT " << xPos << " " << yPos << " " << zPos << std::endl;
        //std::cout << "#END " << xEnd << " " << yEnd << " " << zEnd << std::endl;
        
        //std::cout << xPos << " " << yPos << " " << zPos << std::endl;
        
        // DETERMINE PATH //
        while(Distance>0){
            
            // GET STEP //
            INT Increment=GetStep(xStart,yStart,zStart,xEnd,yEnd,zEnd,xPos,yPos,zPos);
            
            // MAKE STEP //
            if(Increment==-1){
                GoXMinus(xPos,yPos,zPos,1,W,U);  xPos=MOD((xPos-1),Lattice::N[0]);
            }
            if(Increment==+1){
                GoXPlus(xPos,yPos,zPos,1,W,U); xPos=MOD((xPos+1),Lattice::N[0]);
            }
            if(Increment==-2){
                GoYMinus(xPos,yPos,zPos,1,W,U);  yPos=MOD((yPos-1),Lattice::N[1]);
            }
            if(Increment==+2){
                GoYPlus(xPos,yPos,zPos,1,W,U); yPos=MOD((yPos+1),Lattice::N[1]);
            }
            if(Increment==-3){
                GoZMinus(xPos,yPos,zPos,1,W,U);  zPos=MOD((zPos-1),Lattice::N[2]);
            }
            if(Increment==+3){
                GoZPlus(xPos,yPos,zPos,1,W,U); zPos=MOD((zPos+1),Lattice::N[2]);
            }
            
            // CHECK DISTANCE //
            Distance=GetSqrDistance(xPos,yPos,zPos,xEnd,yEnd,zEnd);
            
            // OUTPUT PATH //
            //std::cout << xPos << " " << yPos << " " << zPos << " " << Distance << std::endl;
            
        }
        
        
    }
    
    // DETEMINE TRANSPORTER FOR PRE-DEFINED POSITION AND PATH -- NOTE POSITION GETS UPDATED //
    void GetTransporter(INT &xPos,INT &yPos,INT &zPos,SU_Nc_FUNDAMENTAL_FORMAT *W,GaugeLinks *U,std::vector<INT> Path){
        
        // RESET LINK //
        COPY_SUNcMatrix(W,SUNcGroup::UnitMatrix);
        
        // COMPUTE TRANSPORTER ALONG PATH //
        INT PathLength=Path.size();
        
        for(INT step=0;step<PathLength;step++){
            
            INT Increment=Path[step];
            
            // MAKE STEP //
            if(Increment==-1){
                GoXMinus(xPos,yPos,zPos,1,W,U);  xPos=MOD((xPos-1),Lattice::N[0]);
            }
            if(Increment==+1){
                GoXPlus(xPos,yPos,zPos,1,W,U); xPos=MOD((xPos+1),Lattice::N[0]);
            }
            if(Increment==-2){
                GoYMinus(xPos,yPos,zPos,1,W,U);  yPos=MOD((yPos-1),Lattice::N[1]);
            }
            if(Increment==+2){
                GoYPlus(xPos,yPos,zPos,1,W,U); yPos=MOD((yPos+1),Lattice::N[1]);
            }
            if(Increment==-3){
                GoZMinus(xPos,yPos,zPos,1,W,U);  zPos=MOD((zPos-1),Lattice::N[2]);
            }
            if(Increment==+3){
                GoZPlus(xPos,yPos,zPos,1,W,U); zPos=MOD((zPos+1),Lattice::N[2]);
            }
            
        }

        
    }
    
    // DETEMINE INVERSE TRANSPORTER FOR PRE-DEFINED POSITION AND PATH -- NOTE POSITION GETS UPDATED //
    void GetInverseTransporter(INT &xPos,INT &yPos,INT &zPos,SU_Nc_FUNDAMENTAL_FORMAT *W,GaugeLinks *U,std::vector<INT> Path){
        
        // DETERMINE INVERSE PATH //
        std::vector<INT> InversePath;
        
        GetInversePath(Path,InversePath);
        
        // GET TRANSPORTER ALONG INVERSE PATH //
        GetTransporter(xPos,yPos,zPos,W,U,InversePath);
        
    }

}