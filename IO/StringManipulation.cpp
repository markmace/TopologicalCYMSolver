namespace StringManipulation{
    
    template<typename GenericArgument>
    std::string StringCast(GenericArgument x){
        
        std::stringstream ss;
        ss << x;
        return ss.str();
    }

    template<typename GenericXArgument,typename GenericYArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y){
                       
        std::stringstream ss;
        ss << x << y;
        return ss.str();
        
    }
    
    template<typename GenericXArgument,typename GenericYArgument,typename GenericZArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y,GenericZArgument z){
        
        std::stringstream ss;
        ss << x << y << z;
        return ss.str();
        
    }
    
    template<typename GenericXArgument,typename GenericYArgument,typename GenericZArgument,typename GenericUArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y,GenericZArgument z,GenericUArgument u){
        
        std::stringstream ss;
        ss << x << y << z << u;
        return ss.str();
        
    }
    
    template<typename GenericXArgument,typename GenericYArgument,typename GenericZArgument,typename GenericUArgument,typename GenericVArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y,GenericZArgument z,GenericUArgument u,GenericVArgument v){
        
        std::stringstream ss;
        ss << x << y << z << u << v;
        return ss.str();
        
    }
    
    template<typename GenericXArgument,typename GenericYArgument,typename GenericZArgument,typename GenericUArgument,typename GenericVArgument,typename GenericWArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y,GenericZArgument z,GenericUArgument u,GenericVArgument v,GenericWArgument w){
        
        std::stringstream ss;
        ss << x << y << z << u << v << w;
        return ss.str();
        
    }
    
    template<typename GenericXArgument,typename GenericYArgument,typename GenericZArgument,typename GenericUArgument,typename GenericVArgument,typename GenericWArgument,typename GenericAArgument>
    std::string StringCast(GenericXArgument x,GenericYArgument y,GenericZArgument z,GenericUArgument u,GenericVArgument v,GenericWArgument w,GenericAArgument a){
        
        std::stringstream ss;
        ss << x << y << z << u << v << w << a;
        return ss.str();
        
    }
    
    
}
