        
        //COMPUTE THE DETERMINANT OF A COMPLEX 3x3 MATRIX
        COMPLEX GeneralDet(COMPLEX A[3][3]){
            return A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])+A[0][1]*(A[1][2]*A[2][0] -A[1][0]*A[2][2])+A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
        }
        
        //COMPUTE THE INVERSE OF A COMPLEX 3x3 MATRIX
        void GeneralInverse(COMPLEX A[3][3],COMPLEX AInv[3][3]){
            
            COMPLEX OneOverDetA=DOUBLE(1.0)/GeneralDet(A);
            
            AInv[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])*OneOverDetA;
            AInv[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])*OneOverDetA;
            AInv[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])*OneOverDetA;
            
            AInv[1][0]=(A[1][2]*A[2][0]-A[1][0]*A[2][2])*OneOverDetA;
            AInv[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])*OneOverDetA;
            AInv[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])*OneOverDetA;
            
            AInv[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])*OneOverDetA;
            AInv[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])*OneOverDetA;
            AInv[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])*OneOverDetA;
        }

/*
 namespace Extended{
 
 void Cast(SU_Nc_FUNDAMENTAL_FORMAT *V,COMPLEX VMat[3][3]){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 VMat[i][j]=V[MatrixIndex(i,j)];
 }
 }
 }
 
 void Cast(COMPLEX VMat[3][3],SU_Nc_FUNDAMENTAL_FORMAT *V){
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 V[MatrixIndex(i,j)]=VMat[i][j];
 }
 }            
 }
 
 //COMPUTE THE INVERSE SQUARE ROOT OF A SELF-ADJOINT (HERMITIAN) MATRIX
 void InverseSquareRoot(SU_Nc_FUNDAMENTAL_FORMAT *Q,SU_Nc_FUNDAMENTAL_FORMAT *InvSqrtQ){
 
 //MATRIX BUFFER
 COMPLEX QMat[3][3];
 
 //CAST TO 3x3 MATRIX
 Cast(Q,QMat);
 
 //COMPUTE EIGENSYSTEM OF Q
 COMPLEX U[3][3];
 DOUBLE eV[3];
 
 
 //ONLY WORKS RELIABLY WITH QL FACTORIZATION ALGORITHM
 int SUCCESS=Diagonalization3x3::EigensystemQL(QMat, U, eV);
 
 if(SUCCESS==-1){
 std::cerr << "DIAGONALIZAITON FAILED!" << std::endl;
 exit(0);
 }
 
 //SET BUFFER MATRIX TO ZERO
 for(int i=0;i<3;i++){
 for(int k=0;k<3;k++){
 InvSqrtQ[MatrixIndex(i,k)]=DOUBLE(0.0);
 }
 }
 
 //COMPUTE SQUARE ROOT MATRIX Q=Sqrt(VD.V)
 for(int j=0;j<3;j++){
 
 //COMPUTE PRINCIPAL SQUARE ROOT OF EIGENVALUE
 COMPLEX InvSqrtEval=DSQRT(DOUBLE(1.0)/eV[j]);
 
 //CONTRACT WITH UNITRARY TRANSOFRMATION MATRIX
 for(int i=0;i<3;i++){
 for(int k=0;k<3;k++){
 InvSqrtQ[MatrixIndex(i,k)]+=U[i][j]*InvSqrtEval*conj(U[k][j]);
 }
 }
 
 }
 
 
 }
 
 //COMPUTE UNITARIZATION OF V DEFINED AS VTilde=V.(Sqrt(VD.V))^{-1} * (Det(V)/|Det(V)|)^(-1/3)  //
 void Unitarize(SU_Nc_FUNDAMENTAL_FORMAT *V,SU_Nc_FUNDAMENTAL_FORMAT *VTilde){
 
 SU_Nc_FUNDAMENTAL_FORMAT Q[SUNcGroup::MatrixSize];
 SU_Nc_FUNDAMENTAL_FORMAT K[SUNcGroup::MatrixSize];
 
 //COMPUTE VD.V //
 SUNcGroup::Operations::DU(V,V,Q);	 
 
 //COMPUTE 1/Sqrt(VD.V) //
 InverseSquareRoot(Q,K);
 
 //COMPUTE V.1/Sqrt(VD.V) //
 SUNcGroup::Operations::UU(V,K,VTilde);
 
 //COMPUTE (Det(V)/|Det(V)|)^(-1/3) //
 DOUBLE Arg=-arg( SUNcGroup::Operations::Det(V))/DOUBLE(3.0);
 COMPLEX PhaseFactor=COMPLEX(cos(Arg),sin(Arg));
 
 //COMPUTE V.1/Sqrt(VD.V) (Det(V)/|Det(V)|)^(-1/3) //
 for(int i=0;i<3;i++){
 for(int j=0;j<3;j++){
 VTilde[MatrixIndex(i,j)]=VTilde[MatrixIndex(i,j)]*PhaseFactor;
 }
 }
 
 }       
 }*/
