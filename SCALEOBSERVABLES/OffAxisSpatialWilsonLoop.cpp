#ifndef __OFF_AXIS_SPATIAL_WILSON_LOOP__CPP__
#define __OFF_AXIS_SPATIAL_WILSON_LOOP__CPP__

namespace ScaleObservables {

    namespace OffAxisSpatialWilsonLoop{
        
        DOUBLE CalculateWilsonLoop(INT x1,INT y1,INT z1,INT x2,INT y2,INT z2,INT x3,INT y3,INT z3,GaugeLinks *U){
            
            // BUFFERS FOR WILSON LOOP SIDES //
            SU_Nc_FUNDAMENTAL_FORMAT U12[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT U23[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT U34[SUNcGroup::MatrixSize];
            SU_Nc_FUNDAMENTAL_FORMAT U41[SUNcGroup::MatrixSize];
            
            // DETERMINE PATHES //
            std::vector<INT> p12;   GetPath(x1,y1,z1,x2,y2,z2,p12);
            std::vector<INT> p23;   GetPath(x2,y2,z2,x3,y3,z3,p23);
            
            // GET TRANSPORTERS ALONG THE DIFFERENT PATHES //
            INT xPos=x1; INT yPos=y1; INT zPos=z1;
            
            GetTransporter(xPos,yPos,zPos,U12,GLinks::U,p12);
            GetTransporter(xPos,yPos,zPos,U23,GLinks::U,p23);
            
            GetInverseTransporter(xPos,yPos,zPos,U34,GLinks::U,p12);
            GetInverseTransporter(xPos,yPos,zPos,U41,GLinks::U,p23);
            
            // CHECK THAT LOOP IS CLOSE //
            if(xPos!=x1 || yPos!=y1 || zPos!=z1){
                std::cerr << "#ERROR -- WILSON LOOP IS NOT CLOSED" << std::endl;
                exit(0);
            }
            
            // COMPUTE WILSON LOOP //
            SU_Nc_FUNDAMENTAL_FORMAT ULoop[SUNcGroup::MatrixSize];
            
            SUNcGroup::AdvancedOperations::UUUU(U12,U23,U34,U41,ULoop);
            
            
            // COMPUTE REAL PART OF THE TRACE //
            return SUNcGroup::Operations::ReTr(ULoop)/DOUBLE(Nc);
            
        }

        void WilsonLoopHistogram(std::string fname,GaugeLinks *U){
            
            // COMMANDLINE READOUT
            std::cerr << "#BEGIN WILSON LOOP CALCUALTIONS" << std::endl;
            
            // SETUP HISTOGRAM //
            DOUBLE MAX_LENGTH=std::pow((Lattice::N[0]/2*Lattice::a[0])*(Lattice::N[1]/2*Lattice::a[1])*(Lattice::N[2]/2*Lattice::a[2]),1.0/3.0);
            
            INT NumberOfBins=8*INT(MAX_LENGTH);
            
            Histogram *WLHistogram=new Histogram(0.0,MAX_LENGTH,NumberOfBins);
            
            // MAXIMUM NUMBER OF WILSON LOOPS //
            INT nMax=Lattice::N[0]*Lattice::N[1]*Lattice::N[2];
            
            // COMPUTE nMAX LOOPS //
            for(INT n=0;n<nMax;n++){
                
                ////////////////
                // GET SQUARE //
                ////////////////
                
                // VERTICES OF THE SQUARE //
                INT x1,y1,z1,x2,y2,z2,x3,y3,z3;
                
                // GET SQUARE //
                ScaleObservables::GetSquare(x1,y1,z1,x2,y2,z2,x3,y3,z3,MAX_LENGTH*MAX_LENGTH);
                
                // COMPUTE TRACE //
                DOUBLE WilsonLoopTrace=CalculateWilsonLoop(x1,y1,z1,x2,y2,z2,x3,y3,z3,U);
                
                // CONVERT TO PHYSICAL LENGTH //
                DOUBLE Length=GetDistance(x1,y1,z1,x2,y2,z2)*(U->a[0])*Qs;
                
                // COUNT TO HISTOGRAM //
                WLHistogram->Count(Length,WilsonLoopTrace);
                
            }
            
            // SET OUTPUT FILE//
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            // SET OUTPUT HEADER //
            std::string HeaderMessage=StringManipulation::StringCast("#L ","WL ");
            
            // OUTPUT WILSON LOOP //
            WLHistogram->Output(HeaderMessage,OutputFile);
            
            // CLEAN-UP //
            delete WLHistogram;
            
            // COMMANDLINE READOUT //
            std::cerr << "#FINSIHED WILSON LOOP CALCUALTIONS" << std::endl;
        }
        
        void WilsonLoopHistogram(std::string fname){
            
            WilsonLoopHistogram(fname,GLinks::U);
            
        }
        
        
    }
    
}

#endif
