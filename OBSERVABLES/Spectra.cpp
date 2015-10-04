#ifndef __SPECTRA_CPP__
#define __SPECTRA_CPP__

///////////////////////////////////////////
//INCLUDE FREE FIELD POLARIZATION VECTORS//
///////////////////////////////////////////

#if METRIC_FLAG==BJORKEN_FLAG
#include "../MISC/BjorkenWaveVectors.cpp"
#endif

#if METRIC_FLAG==MINKOWSKI_FLAG
#include "../MISC/MinkowskiWaveVectors.cpp"
#endif

////////////////////////
//SPECTRA MEASUREMENT //
////////////////////////

namespace Observables {
    
    namespace Spectra {
        
        // SET GAUGE FIELDS //
        void SetCoordinateSpaceGaugeFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U){
            
            #pragma omp parallel
            {
                DOUBLE A0Buffer[SUNcAlgebra::VectorSize];
                DOUBLE A1Buffer[SUNcAlgebra::VectorSize];
                DOUBLE A2Buffer[SUNcAlgebra::VectorSize];
                
                #pragma omp for
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            //COMPUTE GAUGE FIELDS
                            SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0)/(U->a[0]),U->Get(x,y,z,0),A0Buffer);
                            SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0)/(U->a[1]),U->Get(x,y,z,1),A1Buffer);
                            SUNcGroup::Operations::ReTrIGenU(DOUBLE(2.0)/(U->a[2]),U->Get(x,y,z,2),A2Buffer);
                            
                            //CAST TO FFTW COMPLEX AND WRITE TO FFTW ARRAY 
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                
                                FourierSpace::A0->SetX(x,y,z,a,A0Buffer[a]);
                                FourierSpace::A1->SetX(x,y,z,a,A1Buffer[a]);
                                FourierSpace::A2->SetX(x,y,z,a,A2Buffer[a]);
                                
                            }
                            
                        }
                    }
                }
                
            } // END PARALLEL
            
        }
        
        // SET ELECTRIC FIELDS //
        void SetCoordinateSpaceElectricFields(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,ElectricFields *E){
                        
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){ 
                        
                        //CAST TO FFTW COMPLEX AND WRITE TO FFTW ARRAY 
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            
                            FourierSpace::E0->SetX(x,y,z,a,(E->a[0]/E->aCube)*(E->Get(x,y,z,0,a)[0]));
                            FourierSpace::E1->SetX(x,y,z,a,(E->a[1]/E->aCube)*(E->Get(x,y,z,1,a)[0]));
                            FourierSpace::E2->SetX(x,y,z,a,(E->a[2]/E->aCube)*(E->Get(x,y,z,2,a)[0]));
                            
                        }
                        
                    }
                } // END PARALLEL FOR
            }
            
        }
        
        
        
        // PERFORM FOURIER TRANSFORM OF GAUGE FIELDS AND ELECTRIC FIELDS //
        void PerformFourierTransform(){
            
            FourierSpace::A0->ExecuteXtoP();  FourierSpace::A1->ExecuteXtoP();  FourierSpace::A2->ExecuteXtoP();
            FourierSpace::E0->ExecuteXtoP();  FourierSpace::E1->ExecuteXtoP();  FourierSpace::E2->ExecuteXtoP();
            
        }
        
        
        // COMPUTE SINGLE PARTICLE OCCUPANCY //
        DOUBLE GetOccupancy(INT pXIndex,INT pYIndex,INT pZIndex, GaugeLinks *U){
            
            // DETERMINE LATTICE MOMENTA //
            COMPLEX cDpx,cDpy,cDpz;
            
            GetMomenta(pXIndex,pYIndex,pZIndex,cDpx,cDpy,cDpz);
            
            // DETERMINE POLARIZATION VECTORS //
            COMPLEX Xi1[Lattice::Dimension]; COMPLEX Xi1Dot[Lattice::Dimension];
            COMPLEX Xi2[Lattice::Dimension]; COMPLEX Xi2Dot[Lattice::Dimension];
            
            PolarizationVectors::Compute(pXIndex,cDpx,pYIndex,cDpy,pZIndex,cDpz,Xi1,Xi1Dot,Xi2,Xi2Dot);
            
            // MODE OCCUPANCY //
            DOUBLE NOcc=0.0;
            
            // GAUGE FIELD MOMENTUM MODES //
            COMPLEX Ap[Lattice::Dimension]; COMPLEX ApDot[Lattice::Dimension];
            
            // COMPUTE OCCUPANCY OF EACH COLOR MODE //
            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                
                // DETERMINE MOMENTUM SPACE GAUGE FIELD //
                Ap[0]=U->aCube*(FourierSpace::A0->GetP(pXIndex,pYIndex,pZIndex,a));
                Ap[1]=U->aCube*(FourierSpace::A1->GetP(pXIndex,pYIndex,pZIndex,a));
                Ap[2]=U->aCube*(FourierSpace::A2->GetP(pXIndex,pYIndex,pZIndex,a));
                
                // DETERMINE TIME DERIVATIVE OF THE MOMENTUM SPACE GAUGE FIELD //
                ApDot[0]=U->aCube*(FourierSpace::E0->GetP(pXIndex,pYIndex,pZIndex,a))*(Dynamics::gDownMetric[0]/Dynamics::MetricDeterminant);
                ApDot[1]=U->aCube*(FourierSpace::E1->GetP(pXIndex,pYIndex,pZIndex,a))*(Dynamics::gDownMetric[1]/Dynamics::MetricDeterminant);
                ApDot[2]=U->aCube*(FourierSpace::E2->GetP(pXIndex,pYIndex,pZIndex,a))*(Dynamics::gDownMetric[2]/Dynamics::MetricDeterminant);
                
                
                // COMPUTE OCCUPATION NUMBER ACCORDING TO PROJECTION FORMULA //
                COMPLEX Xi1Projection=0.0; COMPLEX Xi2Projection=0.0;
                
                for(INT mu=0;mu<Lattice::Dimension;mu++){
                    
                    Xi1Projection+=Dynamics::gUpMetric[mu]*(Xi1Dot[mu]*conj(Ap[mu])-Xi1[mu]*conj(ApDot[mu]));
                    Xi2Projection+=Dynamics::gUpMetric[mu]*(Xi2Dot[mu]*conj(Ap[mu])-Xi2[mu]*conj(ApDot[mu]));
                }
                
                // UPDATE OCCUPATION NUMBER //
                NOcc+=SQR_ABS(Xi1Projection)+SQR_ABS(Xi2Projection);
                
            }
            
            return SQR(Dynamics::MetricDeterminant)*NOcc/(U->Volume*U->aCube);
        }
        
        
        
        // COMPUTE SPECTRA //
        void Compute(std::string fname,GaugeLinks *U,ElectricFields *E){
            
            std::cerr << "#MEASURING SINGLE PARTICLE SPECTRA" << std::endl;
            
            /////////////////////////////////////////////////////
            //      COMPUTE GAUGE FIELDS AND ELECTRIC FIELDS   //
            /////////////////////////////////////////////////////
            
            SetCoordinateSpaceGaugeFields(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U);
            SetCoordinateSpaceElectricFields(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,E);
            
            /////////////////////////////////////////////////////
            //          PERFORM THE FOURIER TRANSFORMS         //
            /////////////////////////////////////////////////////
            
            PerformFourierTransform();
            
            /////////////////////////////////////////////////////
            // COMPUTE THE OBSERVABLE AND WRITE TO OUTPUT FILE //
            /////////////////////////////////////////////////////
            
            // OPEN OUTPOUT STREAM //
            
            Histogram *GluonSpectrum=new Histogram(0.0,2.0*D_SQRT3/Lattice::a[0],8*Lattice::N[0]);
            
            // OUTPUT VARIABLES //
            
            COMPLEX cDpx,cDpy,cDpz; DOUBLE pXValue,pYValue,pZValue; DOUBLE pAbs; DOUBLE ModeOccupancy;
            
            // COMPUTE OCCUPATION NUMBER FOR ALL MOMENTUM MODES //
            for(INT pZIndex=0;pZIndex<=U->N[2]-1;pZIndex++){
                for(INT pYIndex=0;pYIndex<=U->N[1]-1;pYIndex++){
                    for(INT pXIndex=0;pXIndex<=U->N[0]-1;pXIndex++){
                        
                        GetMomenta(pXIndex,pYIndex,pZIndex,cDpx,cDpy,cDpz);
                        
                        pXValue=SIGN(real(cDpx))*sqrt(Dynamics::gUpMetric[0]*SQR(Lattice::aScale)*SQR_ABS(cDpx));
                        pYValue=SIGN(real(cDpy))*sqrt(Dynamics::gUpMetric[1]*SQR(Lattice::aScale)*SQR_ABS(cDpy));
                        pZValue=SIGN(real(cDpz))*sqrt(Dynamics::gUpMetric[2]*SQR(Lattice::aScale)*SQR_ABS(cDpz));
                        
                        ModeOccupancy=GetOccupancy(pXIndex,pYIndex,pZIndex,U)/(2*(Nc*Nc-1));
                        
                        pAbs=sqrt(SQR(pXValue)+SQR(pYValue)+SQR(pZValue));
                        
                        GluonSpectrum->Count(pAbs,ModeOccupancy);
                        
                    }
                    
                    
                }
                
            }
            
            // OUTPUT FILE //
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");

            // CREATE OUPUT //
            GluonSpectrum->Output(OutputFile);
            
            // CLEAN-UP //
            delete GluonSpectrum;
            
        }
     
        void Compute(std::string fname){
            
            Compute(fname,GaugeFixedVariables::U,GaugeFixedVariables::E);
            
        }
        
    }

    
}

#endif
