///////////////////////////////////////////////////////////////////
//   IMPLEMENTATION OF THREE DIMENSIONAL FAST FOURIER TRANSFORM  //
///////////////////////////////////////////////////////////////////

#ifndef __FFT_3D_CPP__
#define __FFT_3D_CPP__

class FFT3D{
    
    // FFTW PARAMETERS //
private:
    
    static const INT MyDimension=3;

    INT MyNx,MyNy,MyNz;

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
    INT XIndex(INT x,INT y,INT z,INT a){
        return a+MyNumberOfTransforms*(x+MyNx*(y+MyNy*z));
    }
    
    INT PIndex(INT px,INT py,INT pz,INT a){
        return a+MyNumberOfTransforms*(px+MyNx*(py+MyNy*pz));
    }

public: 
    
    void SetX(INT x,INT y, INT z,INT a,DOUBLE Value){
        XArray[XIndex(x,y,z,a)][0]=double(Value); XArray[XIndex(x,y,z,a)][1]=double(0.0);
    }
    
    void SetP(INT px,INT py, INT pz,INT a,COMPLEX Value){
        PArray[PIndex(px,py,pz,a)][0]=double(real(Value)); PArray[PIndex(px,py,pz,a)][1]=double(imag(Value));
    }
    
    COMPLEX GetX(INT x,INT y,INT z,INT a){
        return COMPLEX(XArray[XIndex(x,y,z,a)][0],XArray[XIndex(x,y,z,a)][1]);
    }
    
    COMPLEX GetP(INT px,INT py,INT pz,INT a){
        return COMPLEX(PArray[PIndex(px,py,pz,a)][0],PArray[PIndex(px,py,pz,a)][1]);
    }
    
    // FFTW PLANS //
private:
    fftw_plan XtoP;
    fftw_plan PtoX;
    

    
    ///////////////////////
    // CLASS CONSTRUCTOR //
    ///////////////////////
    
public:
    FFT3D(INT Nx,INT Ny,INT Nz,INT NumberOfTransforms){
        
        ////////////////////////
        //SET FFTW PARAMETERS //
        ////////////////////////
        
        //NUMBER OF SITES IN EACH DIMENSION
        MyNx=Nx; MyNy=Ny; MyNz=Nz;
        
        // SET IN REVERSED ORDERING  TO SWITCH FROM ROW-MAJOR TO COLUMN-MAJOR FORMAT //  
        MyNumberOfSites=new INT[3]; 
        
        MyNumberOfSites[0]=MyNz;    MyNumberOfSites[1]=MyNy;    MyNumberOfSites[2]=MyNx;
    
                    
        //NUMBER OF INDEPENDENT TRANSFORMS
        MyNumberOfTransforms=NumberOfTransforms;
        
        //OVERALL DATA SIZE
        MyOverallSize=MyNx*MyNy*MyNz*NumberOfTransforms;
        
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
    
    ~FFT3D(){
        
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
        
        for(INT z=0;z<MyNz;z++){
            for(INT y=0;y<MyNy;y++){
                for(INT x=0;x<MyNx;x++){
                    
                    
                    OutStream << x << " " << y << " " << z << " ";
                    
                    for(INT a=0;a<MyNumberOfTransforms;a++){
                        
                        OutStream << real(GetX(x,y,z,a)) << " " << imag(GetX(x,y,z,a)) << " ";
                    }
                    
                    OutStream << std::endl;
                    
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
        
            for(INT pz=0;pz<MyNz;pz++){
                for(INT py=0;py<MyNy;py++){
                    for(INT px=0;px<MyNx;px++){

                    
                    OutStream << px << " " << py << " " << pz << " ";
                    
                    for(INT a=0;a<MyNumberOfTransforms;a++){
                        
                        OutStream << real(GetP(px,py,pz,a)) << " " << imag(GetP(px,py,pz,a)) << " ";
                    }
                    
                    OutStream << std::endl;
                    
                }
                
                OutStream << std::endl;
                
            }
            
            OutStream << std::endl;
            
        }
        
        OutStream.close();
        
    }
    
};

#endif
