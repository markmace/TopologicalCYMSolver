namespace  IO {
    
    std::string OutputDirectory;
    
    template<typename GenericArgument>
    
    void SetOutputDirectory(GenericArgument x){
        
        OutputDirectory=StringManipulation::StringCast(x,"/");
        
        std::cerr << "#OUTPUT DIRECTORY IS " << OutputDirectory <<  std::endl;
    }

}
