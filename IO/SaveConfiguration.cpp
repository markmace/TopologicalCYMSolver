namespace IO{
    
    void SaveConfiguration(std::string Ufname,std::string Efname,GaugeLinks *U,ElectricFields *E){
        
        // OUTPUT STREAMS //
        std::ofstream UOutStream,EOutStream;
        
        // OUTPUT FILES //
        std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
        std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Efname,"ID",RandomNumberGenerator::MySEED,".txt");

        // OPEN FILES //
        UOutStream.open(UOutputFile.c_str()); EOutStream.open(EOutputFile.c_str());
        
        // SET PRECISION //
        UOutStream.precision(OUTPUT_PRECISION); EOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT z=0;z<U->N[2];z++){
            for(INT y=0;y<U->N[1];y++){
                for(INT x=0;x<U->N[0];x++){
                    for(INT mu=0;mu<Lattice::Dimension;mu++){
                        
                        // OUTPUT GAUGE LINKS //
                        UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(U->Get(x,y,z,mu)) << std::endl;
                        /*
                        // BEGIN SAYANTAN OUTPUT //
                        COMPLEX UMat[Nc*Nc];
                        SUNcGroup::Operations::GetMatrix(GLinks::U->Get(x,y,z,mu),UMat);
                        
                        UOutStream << std::real(UMat[0]) <<  " " << std::imag(UMat[0]) << std::endl;
                        UOutStream << std::real(UMat[2]) <<  " " << std::imag(UMat[2]) << std::endl;
                        UOutStream << std::real(UMat[1]) <<  " " << std::imag(UMat[1]) << std::endl;
                        UOutStream << std::real(UMat[3]) <<  " " << std::imag(UMat[3]) << std::endl;
                        // END SAYANTAN OUTPUT //
                        */
                        // OUTPUT ELECTRIC FIELDS //
                        EOutStream << x << " " << y << " " << z << " " << mu;
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            EOutStream << " " << E->Get(x,y,z,mu,a)[0];
                        }
                        
                        EOutStream << std::endl;
                        
                        
                    }
                    
                }
            }
        }
        
        // CLOSE OUTPUT STREAM //
        UOutStream.close(); EOutStream.close();
        
    }
    // OVERLOAD //
    void SaveConfiguration(std::string Ufname,std::string Efname){
        
        SaveConfiguration(Ufname,Efname,GLinks::U,EFields::E);
    
    }
    
    
}
