#ifndef __MULTI_HISTROGRAM_CPP__
#define __MULTI_HISTROGRAM_CPP__

class MultiHistogram{
    
    DOUBLE *xValues;
    DOUBLE *yValues;
    
    INT NumberOfObservables;
    
    INT *Counts;
    
    INT TotalCounts;
    INT OutOfRange;
    
    INT NBins;
    
    DOUBLE xLow,xHigh;
    
    DOUBLE DeltaX;
    
    public:
    
    void Reset(){
        
        for(INT i=0;i<NBins;i++){
            
            xValues[i]=0.0; Counts[i]=0;
            
            for(INT j=0;j<NumberOfObservables;j++){
                yValues[NumberOfObservables*i+j]=0.0;
            }
            
        }
        
        TotalCounts=0; OutOfRange=0;
        
    }
    
    
    void Count(DOUBLE x,DOUBLE y[]){
        
        INT pos=INT((x-xLow)/DeltaX);
        
        if(pos>=0 && pos<NBins){
            
            xValues[pos]+=x;
            
            for(INT j=0;j<NumberOfObservables;j++){
                yValues[NumberOfObservables*pos+j]+=y[j];
            }
            
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
                
                OutStream << xValues[i]/Counts[i];
                
                for(INT j=0;j<NumberOfObservables;j++){
                    OutStream << " " << yValues[NumberOfObservables*i+j]/Counts[i];
                }
                
                OutStream << " " << Counts[i]/DOUBLE(TotalCounts+OutOfRange) << std::endl;
                
            }
            
        }
        
        OutStream.close();
        
    }
    
    
    MultiHistogram(DOUBLE xL,DOUBLE xH,INT N,INT NObs){
        
        xLow=xL; xHigh=xH; NBins=N; DeltaX=(xHigh-xLow)/NBins;
        
        xValues=new DOUBLE[NBins];  Counts=new INT[NBins];
        
        NumberOfObservables=NObs;
        
        yValues=new DOUBLE[NumberOfObservables*NBins];
        
        Reset();
        
    }

    MultiHistogram(){
      
        delete[] xValues; 
        delete[] yValues;
        delete[] Counts;

        
    }
    
    
    
    
};

#endif
