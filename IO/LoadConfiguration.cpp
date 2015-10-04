namespace IO{
    
    void LoadConfiguration(GaugeLinks *U,ElectricFields *E){
        
        // INPUT STREAMS //
        std::ifstream UInStream,EInStream;
        
        // INPUT FILES //
        std::string UInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputUFile);
        std::string EInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputEFile);
        
        // OPEN FILES //
        UInStream.open(UInputFile.c_str()); EInStream.open(EInputFile.c_str());
        
        // SET PRECISION //
        UInStream.precision(OUTPUT_PRECISION); EInStream.precision(OUTPUT_PRECISION);
        
        // INPUT DATA //
        INT x; INT y; INT z; INT mu;
        DOUBLE U0; DOUBLE U1; DOUBLE U2; DOUBLE U3;
        DOUBLE E0; DOUBLE E1; DOUBLE E2;
        
        // BUFFER FOR STRING TO MATRIX //
        std::string UStr;
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        
        // INPUT GAUGE LINKS //
        if (UInStream.is_open()){
            while(!UInStream.eof()){
                
                // INPUT GAUGE LINKS //
                UInStream >> x >> y >> z >> mu >> U0 >> U1 >> U2 >> U3;
                
                //STRING BUFFER FOR MATRIX
                UStr = StringManipulation::StringCast(U0," ",U1," ",U2," ",U3);
                
                // SET GAUGE LINKS AS MATRIX //
                SUNcGroup::IO::StringToMatrix(UStr,U->Get(x,y,z,mu));
                
            }
            
        }
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // INPUT ELECTRIC FIELDS //
        if (EInStream.is_open()){
            while(!EInStream.eof()){
                
                // INPUT ELECTRIC FIELDS //
                EInStream >> x >> y >> z >> mu >> E0 >> E1 >> E2;
                
                E->Get(x,y,z,mu,0)[0]=E0;
                E->Get(x,y,z,mu,1)[0]=E1;
                E->Get(x,y,z,mu,2)[0]=E2;
                
            }
            
        }
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // CLOSE INPUT STREAM //
        UInStream.close(); EInStream.close();
        
    }
    
    
    void LoadLinks(GaugeLinks *U){
        
        // INPUT STREAMS //
        std::ifstream UInStream;
        
        // INPUT FILES //
        std::string UInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputUFile);
        
        // OPEN FILES //
        UInStream.open(UInputFile.c_str());
        
        // SET PRECISION //
        UInStream.precision(OUTPUT_PRECISION);
        
        // INPUT DATA //
        INT x; INT y; INT z; INT mu;
        DOUBLE U0; DOUBLE U1; DOUBLE U2; DOUBLE U3;
        
        // BUFFER FOR STRING TO MATRIX //
        std::string UStr;
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        
        // INPUT GAUGE LINKS //
        if (UInStream.is_open()){
            while(!UInStream.eof()){
                
                // INPUT GAUGE LINKS //
                UInStream >> x >> y >> z >> mu >> U0 >> U1 >> U2 >> U3;
                
                
                //STRING BUFFER FOR MATRIX
                UStr = StringManipulation::StringCast(U0," ",U1," ",U2," ",U3);
                
                // SET GAUGE LINKS AS MATRIX //
                SUNcGroup::IO::StringToMatrix(UStr,U->Get(x,y,z,mu));
                
            }
            
        }
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        // CLOSE INPUT STREAM //
        UInStream.close();
        
    }
    
    void LoadFields(ElectricFields *E){
        
        // INPUT STREAMS //
        std::ifstream EInStream;
        
        // INPUT FILES //
        std::string EInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputEFile);
        
        // OPEN FILES //
         EInStream.open(EInputFile.c_str());
        
        // SET PRECISION //
        EInStream.precision(OUTPUT_PRECISION);
        
        // INPUT DATA //
        INT x; INT y; INT z; INT mu;
        DOUBLE E0; DOUBLE E1; DOUBLE E2;
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // INPUT ELECTRIC FIELDS //
        if (EInStream.is_open()){
            while(!EInStream.eof()){
                
                // INPUT ELECTRIC FIELDS //
                EInStream >> x >> y >> z >> mu >> E0 >> E1 >> E2;
                
                E->Get(x,y,z,mu,0)[0]=E0;
                E->Get(x,y,z,mu,1)[0]=E1;
                E->Get(x,y,z,mu,2)[0]=E2;
                
            }
            
        }
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // CLOSE INPUT STREAM //
        EInStream.close();
        
    }

    
}
/*
// LOAD GAUGE LINKS AND ELECTRIC FIELDS FROM FILE //
void Load(std::string Ufname,std::string Efname,GaugeLinks *U,ElectricFields *E){
    
    std::cerr << "#LOADING GAUGE LINKS " << Ufname.c_str() << std::endl;
    
    //LOAD FROM FILE
    std::ifstream UIn;
    std::ifstream EIn;
    
    UIn.open(Ufname.c_str());
    EIn.open(Efname.c_str());
    
    std::string ULink;
    std::string EField;

    
    //GET POSITION VALUES
    INT xIN,yIN,zIN,muIN;
    
    //NUMBER OF LINES READ FROM FILE
    INT InputCount=0;
    
    // MONITOR PRECISION //
    DOUBLE MaxUnitarityViolation=DOUBLE(0.0);
    
    //GET GAUGE LINK DATA FROM INPUT FILES
    while(UIn.good()){
        
        //READ LINES
        getline(UIn,ULink);
        
        //PROCESS FILE LINE BY LINE
        if(!(ULink.empty())){
            
            //STRING TOKEN
            std::stringstream ULinkValues(ULink);
            
            //GET POSITIONS IN FILE
            ULinkValues >> xIN; ULinkValues >> yIN; ULinkValues >> zIN; ULinkValues >> muIN;
            
            //CHECK POSITIONS AND SET VALUES TO GAUGE LINK ARRAY
            if(InputCount==U->Index(xIN,yIN,zIN,muIN)){
                
                std::stringstream UString;    std::string strBuff;
                
                while(ULinkValues >> strBuff){UString << strBuff << " ";}
                
                SUNcGroup::IO::StringToMatrix(UString.str(),U->Get(xIN,yIN,zIN,muIN));
                
                MaxUnitarityViolation=std::max(MaxUnitarityViolation,fabs(SUNcGroup::Operations::UnitarityNorm(U->Get(xIN,yIN,zIN,muIN))));
                
            }
            
            //BREAK IF POSITIONS DO NOT MATCH
            else{
                std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- TRANSVERSE POSITIONS DO NOT MATCH" << std::endl;
                exit(0);
            }
            
            //INCREASE POSITION COUNT
            InputCount++;
        }
    }
    
    std::cerr << "#LOADING ELECTRIC FIELDS " << Efname.c_str() << std::endl;

    //GET ELECTRIC FIELD DATA FROM INPUT FILES
    while(EIn.good()){
        
        //READ LINES
        getline(EIn,EField);
        
        //PROCESS FILE LINE BY LINE
        if(!(EField.empty())){
            
            //STRING TOKENIZE
            std::stringstream EFieldValues(EField);
            
            //GET POSITIONS IN FILE
            EFieldValues >> xIN; EFieldValues >> yIN; EFieldValues >> zIN;EFieldValues >> muIN;EFieldValues >> aIN;
            
            //CHECK POSITIONS AND SET VALUES TO ELECTRIC FIELD ARRAY
            if(InputCount==E->Index(xIN,yIN,zIN,muIN,aIN)){
                
                std::stringstream UString;    std::string strBuff;
                
                while(EFieldValues >> strBuff){EString << strBuff << " ";}
                
                E->Get(x,y,z,mu,0)[0]=E0;
                E->Get(x,y,z,mu,1)[0]=E1;
                E->Get(x,y,z,mu,2)[0]=E2;
 
            }
 
            //BREAK IF POSITIONS DO NOT MATCH
            else{
                std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- TRANSVERSE POSITIONS DO NOT MATCH" << std::endl;
                exit(0);
            }
            
            //INCREASE POSITION COUNT
            InputCount++;
        }
    }
    
    if(InputCount!=Lattice::N[0]*Lattice::N[1]*Lattice::N[2]){
        std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- GAUGE LINKS NOT LOADED CORRECTLY" << std::endl;
        exit(0);
    }
    
    std::cerr << "#MAX UNITARITY VIOLATION IS " << MaxUnitarityViolation << std::endl;
    
}
*/

