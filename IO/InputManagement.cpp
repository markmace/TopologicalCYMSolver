namespace  IO {
    
    std::string InputDirectory;
    
    template<typename GenericArgument>
    
    void SetInputDirectory(GenericArgument x){
        
        InputDirectory=StringManipulation::StringCast(x,"/");
        
        std::cerr << "#INPUT DIRECTORY IS " << InputDirectory <<  std::endl;
    }
    
    std::string InputUFile;
    std::string InputEFile;
    
    template<typename GenericArgument>
    
    void SetInputFile(GenericArgument x,GenericArgument y){
        
        InputUFile=StringManipulation::StringCast("UOutT",x,"ID",y,".txt");
        
        std::cerr << "#INPUT GAUGE LINK FILE IS " << InputUFile <<  std::endl;
        
        InputEFile=StringManipulation::StringCast("EOutT",x,"ID",y,".txt");
        
        std::cerr << "#INPUT ELECTRIC FIELD FILE IS " << InputEFile <<  std::endl;
    }
    
}
