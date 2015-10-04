#ifndef __HISTROGRAM_CPP__
#define __HISTROGRAM_CPP__

class Histogram{
    
    DOUBLE *xValues;
    DOUBLE *yValues;
    
    INT *Counts;
    
    INT TotalCounts;
    INT OutOfRange;
    
    INT NBins;
    
    DOUBLE xLow,xHigh;
    
    DOUBLE DeltaX;
    
    public:
    
    void Reset(){
        
        for(INT i=0;i<NBins;i++){
            
            xValues[i]=0.0; yValues[i]=0.0; Counts[i]=0; 
            
        }
        
        TotalCounts=0; OutOfRange=0;
        
    }
    
    
    void Count(DOUBLE x,DOUBLE y){
        
        INT pos=INT((x-xLow)/DeltaX);
        
        if(pos>=0 && pos<NBins){
            
            xValues[pos]+=x;
            yValues[pos]+=y;
            
            Counts[pos]++;
            TotalCounts++;
        }
        
        else{
            OutOfRange++;
        }
    }
    
    void Output(std::string fname){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        OutStream << "#xLow=" << xLow << " xHigh=" << xHigh << " DeltaX=" << DeltaX << " NBins=" << NBins << std::endl;
        OutStream << "#TOTALCOUNTS=" << TotalCounts << " UNCOUNTED=" << OutOfRange << std::endl;
        
        for(INT i=0;i<NBins;i++){
            
            if(Counts[i]>0){
                
                OutStream << xValues[i]/Counts[i] << " " << yValues[i]/Counts[i] << " " << Counts[i]/DOUBLE(TotalCounts+OutOfRange) << std::endl;
                
            }
            
        }
        
        OutStream.close();
        
    }
    
    
    Histogram(DOUBLE xL,DOUBLE xH,INT N){
        
        xLow=xL; xHigh=xH; NBins=N; DeltaX=(xHigh-xLow)/NBins;
        
        xValues=new DOUBLE[NBins]; yValues=new DOUBLE[NBins]; Counts=new INT[NBins];
        
        Reset();
        
    }

    ~Histogram(){
      
        delete[] xValues; 
        delete[] yValues;
        delete[] Counts;
        
    }
    
    
    
    
};

#endif
