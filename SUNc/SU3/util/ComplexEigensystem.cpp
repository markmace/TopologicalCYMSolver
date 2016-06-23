#ifndef __COMPLEX__EIGENSYSTEM__CPP__
#define __COMPLEX__EIGENSYSTEM__CPP__

/* ZGEEV prototype */
extern "C" {
    void zgeev(char const* jobvl, char const* jobvr, int* n, COMPLEX* a, int* lda, COMPLEX* w, COMPLEX* vl, int* ldvl,COMPLEX* vr, int* ldvr, COMPLEX* work, int* lwork, double* rwork, int* info );
}

/* Parameters */
#define LDA Nc
#define LDVL Nc
#define LDVR Nc

////////////////////////////////////////
//                                    //
// NOTE: LAPACK DEFINES MATRICES AS   //
//                                    //
//          A_0 A_3 A_6               //
//          A_1 A_4 A_7               //
//          A_2 A_5 A_8               //
//                                    //
////////////////////////////////////////

namespace ComplexEigensystem{
    
    INT GetComplexEigensystem(COMPLEX A[9],COMPLEX EVals[3], COMPLEX EVecs[9]) {
        
        /* Locals */
        int n = Nc, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
        COMPLEX wkopt;
        COMPLEX* work;
        /* Local arrays */
        /* rwork dimension should be at least 2*n */
        double rwork[2*Nc];
        COMPLEX vl[Nc*Nc]; // LEFT EIGENVECTORS
        
        /* Executable statements */
        /* Query and allocate the optimal workspace */
        lwork = -1;
        zgeev( "Vectors", "Vectors", &n, A, &lda, EVals, vl, &ldvl, EVecs, &ldvr, &wkopt, &lwork, rwork, &info );
        lwork = int(real(wkopt));
        work = (COMPLEX*)malloc( lwork*sizeof(COMPLEX) );
        /* Solve eigenproblem */
        zgeev( "Vectors", "Vectors", &n, A, &lda, EVals, vl, &ldvl, EVecs, &ldvr,work, &lwork, rwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            return -1;
            exit( 1 );
        }
        
        /* Free workspace */
        free( (void*)work );
        
        return 0;
    }
}

#endif

