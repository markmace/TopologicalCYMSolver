namespace Observables{
    
    namespace Correlators{
        
        FFT3D *Correlators;
        
        void Init(){
            
            // ALLOCATE FOURIER TRANSFORM OF E,B,E.B SQUARED CORRELATOR //
            Correlators=new FFT3D(Lattice::N[0],Lattice::N[1],Lattice::N[2],3);
            
        }
        
        void Measure(std::string fname,GaugeLinks *U,ElectricFields *E){
            
            std::cerr << "#MEASURING E AND B CORREALTORS" << std::endl;
            
            // NORMALIZATION FACTOR //
            DOUBLE NormalizationFactor=1.0/(Lattice::N[0]*Lattice::N[1]*Lattice::N[2]);
            
            //BUFFERS FOR AVERAGE FIELD STRENGTH
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            // CONSTANTS NEEDED TO E AND B SQUARED //
            DOUBLE cB[Lattice::Dimension];
            DOUBLE cE[Lattice::Dimension];
            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                cB[mu]=(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                cE[mu]=(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * SQR(U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
            
            // CONSTANTS NEEDED FOR E.B SQUARED //
            DOUBLE cS[Lattice::Dimension];
            
            for(INT mu=0;mu<Lattice::Dimension;mu++){
                cS[mu]=Lattice::aScale*SQR(U->a[mu])/U->aCube;
            }
            
            //SET LOCAL VALUES
            DOUBLE EELoc,BBLoc,EBEBLoc;
            
            for(INT z=0;z<Lattice::N[2];z++){
                for(INT y=0;y<Lattice::N[1];y++){
                    for(INT x=0;x<Lattice::N[0];x++){
                        
                        COMPUTE_AVG_FIELD_STRENGTH(x,y,z);
                        
                        // UPDATE E.B //
                        EELoc=cE[0]*E0SqrLoc+cE[1]*E1SqrLoc+cE[2]*E2SqrLoc;
                        BBLoc=cB[0]*B0SqrLoc+cB[1]*B1SqrLoc+cB[2]*B2SqrLoc;
                        EBEBLoc=cS[0]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B0Loc)+cS[1]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B1Loc)+cS[2]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B2Loc);
                        // SET CORRELATOR IN POSITION SPACE INITALLY //
                        Correlators->SetX(x,y,z,0,EELoc);
                        Correlators->SetX(x,y,z,1,BBLoc);
                        Correlators->SetX(x,y,z,2,EBEBLoc);
                        
                        
                    }
                }
            }

            // PERFORM FFT //
            Correlators->ExecuteXtoP();
            
            // COMPUTE CORRELATORS IN FOURIER SPACE //
            MultiHistogram *MomentumCorrelators=new MultiHistogram(DOUBLE(0.0),2.0*D_SQRT3,4*Lattice::N[0],3);
            
            COMPLEX cDpx,cDpy,cDpz; DOUBLE pAbs;
            
            for(INT pZIndex=0;pZIndex<Lattice::N[2];pZIndex++){
                for(INT pYIndex=0;pYIndex<Lattice::N[1];pYIndex++){
                    for(INT pXIndex=0;pXIndex<Lattice::N[0];pXIndex++){
                        
                        // COMPUTE FOURIER TRACE //
                        COMPLEX E2=COMPLEX(0.0,0.0);
                        COMPLEX B2=COMPLEX(0.0,0.0);
                        COMPLEX EB2=COMPLEX(0.0,0.0);
                        
                        COMPLEX EE=Correlators->GetP(pXIndex,pYIndex,pZIndex,0);
                        COMPLEX BB=Correlators->GetP(pXIndex,pYIndex,pZIndex,1);
                        COMPLEX EBEB=Correlators->GetP(pXIndex,pYIndex,pZIndex,2);
                        
                        E2=NormalizationFactor*(EE*conj(EE));
                        B2=NormalizationFactor*(BB*conj(BB));
                        EB2=NormalizationFactor*(EBEB*conj(EBEB));
                        
                        // WRITE TO ARRAY //
                        Correlators->SetP(pXIndex,pYIndex,pZIndex,0,E2);
                        Correlators->SetP(pXIndex,pYIndex,pZIndex,1,B2);
                        Correlators->SetP(pXIndex,pYIndex,pZIndex,2,EB2);
                        
                        GetAbsMomentum(pXIndex,pYIndex,pZIndex,pAbs)
                        
                        DOUBLE Values[3]={abs(E2),abs(B2),abs(EB2)};
                        
                        MomentumCorrelators->Count(pAbs,Values);

                        
                    }
                }
            }
            
            // OUTPUT MOMENTUM SPACE CORRELATORS //
            std::string MomentumCorrelatorsOutputFile=StringManipulation::StringCast(IO::OutputDirectory,"Momentum",fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            MomentumCorrelators->Output(MomentumCorrelatorsOutputFile);
            
            delete MomentumCorrelators;
            
            // PEFRORM FFT //
            Correlators->ExecutePtoX();
            
            // CREATE HISTOGRAM //
            MultiHistogram *CorrelatorsResult=new MultiHistogram(DOUBLE(0.0),(Lattice::N[0]/2),4*Lattice::N[0],3);
            
            // COUNT ONLY HALF OF THE PHASE-SPACE TO ACCOUNT FOR PERIODICITY //
            for(INT z=0;z<Lattice::N[2]/2;z++){
                for(INT y=0;y<Lattice::N[1]/2;y++){
                    for(INT x=0;x<Lattice::N[0]/2;x++){
                        
                        // GET DISTANCE AND E,B,E.B SQUARED AMPLITUDE //
                        DOUBLE r=sqrt(SQR(x)+SQR(y)+SQR(z));
                        DOUBLE De2=NormalizationFactor*real(Correlators->GetX(x,y,z,0)); DOUBLE Db2=NormalizationFactor*real(Correlators->GetX(x,y,z,1)); DOUBLE Deb2=NormalizationFactor*real(Correlators->GetX(x,y,z,2));
                        
                        DOUBLE Values[3]={De2,Db2,Deb2};
                        
                        // COUNT TO HISTOGRAM //
                        CorrelatorsResult->Count(r,Values);

                    }
                }
            }
            
            // SET OUPUT FILE //
            std::string CorrelatorsOutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            // CREATE OUPUT //
            CorrelatorsResult->Output(CorrelatorsOutputFile);
            
            // CLEAN-UP //
            delete CorrelatorsResult;
            
        }
        
        void Measure(std::string fname){
            
            Measure(fname,GLinks::U,EFields::E);
            
        }
        
    }
    


}

