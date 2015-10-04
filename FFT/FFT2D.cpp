///////////////////////////////////////////////////////////////////
//   IMPLEMENTATION OF THREE DIMENSIONAL FAST FOURIER TRANSFORM  //
///////////////////////////////////////////////////////////////////

#ifndef __FFT_2D_CPP__
#define __FFT_2D_CPP__

class FFT2D{
    
    // FFTW PARAMETERS //
private:
    
    static const INT MyDimension=2;
    
    INT MyNx,MyNy;
    
    INT *MyNumberOfSites;
    
    INT MyNumberOfTransforms;
    
    INT MyOverallSize;
    
    INT MyInputMemDistance;
    INT MyOutputMemDistance;    
    
    INT MyInputDataSpacing;
    INT MyOutputDataSpacing;
    
    INT *MyInEmbed;
    INT *MyOutEmbed;
    
    
    // DATA STORAGE //
private:
    fftw_complex *XArray;
    fftw_complex *PArray;
    
    // INDEXING IN COORDINATE AND MOMENTUM SPACE //
    
private:    
    INT XIndex(INT x,INT y,INT a){
        return a+MyNumberOfTransforms*(x+MyNx*y);
    }
    
    INT PIndex(INT px,INT py,INT a){
        return a+MyNumberOfTransforms*(px+MyNx*py);
    }
    
public: 
    
    void SetX(INT x,INT y,INT a,DOUBLE Value){
        XArray[XIndex(x,y,a)][0]=double(Value); XArray[XIndex(x,y,a)][1]=double(0.0);
    }
    
    void SetP(INT px,INT py,INT a,COMPLEX Value){
        PArray[PIndex(px,py,a)][0]=double(real(Value)); PArray[PIndex(px,py,a)][1]=double(imag(Value));
    }
    
    COMPLEX GetX(INT x,INT y,INT a){
        return COMPLEX(XArray[XIndex(x,y,a)][0],XArray[XIndex(x,y,a)][1]);
    }
    
    COMPLEX GetP(INT px,INT py,INT a){
        return COMPLEX(PArray[PIndex(px,py,a)][0],PArray[PIndex(px,py,a)][1]);
    }
    
    // FFTW PLANS //
private:
    fftw_plan XtoP;
    fftw_plan PtoX;
    
    
    
    ///////////////////////
    // CLASS CONSTRUCTOR //
    ///////////////////////
    
public:
    FFT2D(INT Nx,INT Ny,INT NumberOfTransforms){
        
        ////////////////////////
        //SET FFTW PARAMETERS //
        ////////////////////////
        
        //NUMBER OF SITES IN EACH DIMENSION
        MyNx=Nx; MyNy=Ny;
        
        // SET IN REVERSED ORDERING  TO SWITCH FROM ROW-MAJOR TO COLUMN-MAJOR FORMAT //  
        MyNumberOfSites=new INT[2]; 
        
        MyNumberOfSites[0]=MyNy;    MyNumberOfSites[1]=MyNx;
        
        
        //NUMBER OF INDEPENDENT TRANSFORMS
        MyNumberOfTransforms=NumberOfTransforms;
        
        //OVERALL DATA SIZE
        MyOverallSize=MyNx*MyNy*NumberOfTransforms;
        
        //SPACING BETWEEN DATASETS
        MyInputMemDistance=1; MyOutputMemDistance=1;
        
        //SPACING BETWEEN POINTS OF THE SAME DATASET
        MyInputDataSpacing=MyNumberOfTransforms;  MyOutputDataSpacing=MyNumberOfTransforms;
        
        // EMBEDDING OF INDIVIDUAL DATASETS //
        MyInEmbed=MyNumberOfSites;
        MyOutEmbed=MyNumberOfSites;
        
        ////////////////////////
        //  ALLOCATE ARRAYS   //
        ////////////////////////
        
        XArray=fftw_alloc_complex(MyOverallSize);
        PArray=fftw_alloc_complex(MyOverallSize);
        
        ////////////////////////
        //  COMPUTE PLANS     //
        ////////////////////////
        
        
        XtoP=fftw_plan_many_dft(MyDimension,MyNumberOfSites,MyNumberOfTransforms,
                                XArray,MyInEmbed,
                                MyInputDataSpacing, MyInputMemDistance,
                                PArray,MyOutEmbed,
                                MyOutputDataSpacing,MyOutputMemDistance,
                                FFTW_FORWARD,MY_FFTW_PLANNER_FLAG);
        
        PtoX=fftw_plan_many_dft(MyDimension,MyNumberOfSites,MyNumberOfTransforms,
                                PArray,MyOutEmbed,
                                MyOutputDataSpacing, MyOutputMemDistance,
                                XArray,MyInEmbed,
                                MyInputDataSpacing, MyInputMemDistance,
                                FFTW_BACKWARD,MY_FFTW_PLANNER_FLAG);
        
        
    }
    
    ///////////////////////
    // CLASS DESTRUCTOR  //
    ///////////////////////
    
    ~FFT2D(){
        
        //DESTROY PLANS//
        fftw_destroy_plan(XtoP);
        fftw_destroy_plan(PtoX);
        
        //DELETE FFTW PARAMETERS //
        delete[] MyNumberOfSites;
        
        //FREE MEMORY //
        fftw_free(XArray);
        fftw_free(PArray);
        
    }
    
    
    //////////////////////////
    //  EXECUTION COMMANDS  //
    //////////////////////////
public:
    
    void ExecuteXtoP(){
        fftw_execute(XtoP);
    }
    
    void ExecutePtoX(){
        fftw_execute(PtoX);
    }
    
public:
    
    void OutputX(std::string fname){
        
        std::ofstream OutStream;
        
        OutStream.open(fname.c_str());
        
        for(INT y=0;y<MyNy;y++){
            for(INT x=0;x<MyNx;x++){
                
                
                OutStream << x << " " << y << " ";
                
                for(INT a=0;a<MyNumberOfTransforms;a++){
                    
                    OutStream << real(GetX(x,y,a)) << " " << imag(GetX(x,y,a)) << " ";
                }
                
                OutStream << std::endl;
                
            }
            
            OutStream << std::endl;
            
        }
        
        
        OutStream.close();
        
    }
    
    void OutputP(std::string fname){
        
        std::ofstream OutStream;
        
        OutStream.open(fname.c_str());
        
        for(INT py=0;py<MyNy;py++){
            for(INT px=0;px<MyNx;px++){
                
                
                OutStream << px << " " << py << " ";
                
                for(INT a=0;a<MyNumberOfTransforms;a++){
                    
                    OutStream << real(GetP(px,py,a)) << " " << imag(GetP(px,py,a)) << " ";
                }
                
                OutStream << std::endl;
                
            }
            
            OutStream << std::endl;
            
        }
        
        OutStream.close();
        
    }
    
};

#endif
