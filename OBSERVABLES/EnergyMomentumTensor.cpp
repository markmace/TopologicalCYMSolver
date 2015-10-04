#ifndef __ENERGYDENSITYMAP__CPP__
#define __ENERGYDENSITYMAP__CPP__

namespace Observables {
        
    namespace EnergyMomentumTensor{
        
        void CreateMap(std::string fname,INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
            
            //CONSTANTS NEEDED TO COMPUTE ENERGY MOMENTUM TENSOR
            DOUBLE qB[Lattice::Dimension];  DOUBLE qE[Lattice::Dimension];
            
            //SET CONSTANTS
            for(int mu=0;mu<Lattice::Dimension;mu++){
                
                qB[mu]=sqrt(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * (U->a[mu]*SQR(Lattice::aScale)/U->aCube);
                qE[mu]=sqrt(Dynamics::gDownMetric[mu]/SQR(Dynamics::MetricDeterminant)) * (U->a[mu]*SQR(Lattice::aScale)/U->aCube);
            }
            
            //DIAGONAL COMPONENTS OF THE ENERGY MOMENTUM TENSOR
            DOUBLE T00,TXX,TYY,TZZ;
            DOUBLE T0X,T0Y,T0Z,TXY,TYZ,TZX;
                        
            //BUFFERS FOR AVERAGE FIELD STRENGTH
            SET_AVG_FIELD_STRENGTH_BUFFERS();
            
            //CREATE OUTPUT STREAM
            std::ofstream OutStream;
            
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".txt");
            
            OutStream.open(OutputFile.c_str());
            
            //UPDATE ALL GAUGE LINKS
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){

                        //COMPUTE AVERAGE FIELD STRENGTH
                        COMPUTE_AVG_FIELD_STRENGTH(x,y,z);
                        
                        //COMPUTE DIAGONAL COMPONENTS OF ENERGY MOMENTUM TENSOR
                        T00=DOUBLE(0.5)*(SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc);
                        TXX=DOUBLE(0.5)*((SQR(qE[1])*E1SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[1])*B1SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qE[0])*E0SqrLoc+SQR(qB[0])*B0SqrLoc));
                        TYY=DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[2])*E2SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[2])*B2SqrLoc)-(SQR(qE[1])*E1SqrLoc+SQR(qB[1])*B1SqrLoc));
                        TZZ=DOUBLE(0.5)*((SQR(qE[0])*E0SqrLoc+SQR(qE[1])*E1SqrLoc+SQR(qB[0])*B0SqrLoc+SQR(qB[1])*B1SqrLoc)-(SQR(qE[2])*E2SqrLoc+SQR(qB[2])*B2SqrLoc));
                        
                        //COMPUTE POYNTING VECTOR
                        T0X=qE[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B2Loc)-qE[2]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B1Loc);
                        T0Y=qE[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,B0Loc)-qE[0]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B2Loc);
                        T0Z=qE[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,B1Loc)-qE[1]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,B0Loc);

                        //COMPUTE SPATIAL STRESS
                        TXY=qE[0]*qE[1]*SUNcAlgebra::Operations::ScalarProduct(E0Loc,E1Loc)+qB[0]*qB[1]*SUNcAlgebra::Operations::ScalarProduct(B0Loc,B1Loc);
                        TYZ=qE[1]*qE[2]*SUNcAlgebra::Operations::ScalarProduct(E1Loc,E2Loc)+qB[1]*qB[2]*SUNcAlgebra::Operations::ScalarProduct(B1Loc,B2Loc);
                        TZX=qE[2]*qE[0]*SUNcAlgebra::Operations::ScalarProduct(E2Loc,E0Loc)+qB[2]*qB[0]*SUNcAlgebra::Operations::ScalarProduct(B2Loc,B0Loc);
                        
                                                                        
                        //WRITE OUTPUT TO FILE
                        OutStream << x << " " << y << " " << z << " ";
                        
                        OutStream << T00 << " " << TXX << " " << TYY << " " << TZZ << " ";
                        OutStream << T0X << " " << T0Y << " " << T0Z << " ";
                        OutStream << TXY << " " << TYZ << " " << TZX << " ";
                        
                        OutStream << std::endl;
                        
                    }
                    
                    OutStream << std::endl;
                    
                }
                                
            }
            
            //CLOSE OUTPUT STREAM
            OutStream.close();
            
        }
        
        
        //CREATE ENERGY DENSITY MAP
        void CreateMap(std::string fname,GaugeLinks *U,ElectricFields *E){
            
            CreateMap(fname,0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
            
        }
        
        
        void CreateMap(std::string fname){
            
            CreateMap(fname,GLinks::U,EFields::E);
            
        }
        
    }
    
}


#endif
