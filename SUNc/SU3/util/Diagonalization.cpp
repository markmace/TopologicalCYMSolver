#ifndef __DIAGONALIZATION__3x3__CPP__
#define __DIAGONALIZATION__3x3__CPP__

// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------

#include <cfloat>

// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a hermitian 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------

namespace Diagonalization3x3{
    
    int Eigenvalues(COMPLEX A[3][3], DOUBLE w[3]){
        
        DOUBLE m,c1, c0;
        
        // Determine coefficients of characteristic poynomial
        // ----------------------------------------------------------------------------
        //       | a   d   f  |
        //  A =  | d*  b   e  |
        //       | f*  e*  c  |
        // ----------------------------------------------------------------------------
        
        COMPLEX de = A[0][1] * A[1][2];                   // d * e
        
        DOUBLE dd = SQR_ABS(A[0][1]);                                  // d * conj(d)
        DOUBLE ee = SQR_ABS(A[1][2]);                                  // e * conj(e)
        DOUBLE ff = SQR_ABS(A[0][2]);                                  // f * conj(f)
        m  = real(A[0][0]) + real(A[1][1]) + real(A[2][2]);
        c1 = (real(A[0][0])*real(A[1][1])  // a*b + a*c + b*c - d*conj(d) - e*conj(e) - f*conj(f)
              + real(A[0][0])*real(A[2][2])
              + real(A[1][1])*real(A[2][2]))
        - (dd + ee + ff);
        c0 = real(A[2][2])*dd + real(A[0][0])*ee + real(A[1][1])*ff
        - real(A[0][0])*real(A[1][1])*real(A[2][2])
        - DOUBLE(2.0) * (real(A[0][2])*real(de) + imag(A[0][2])*imag(de));
        // c*d*conj(d) + a*e*conj(e) + b*f*conj(f) - a*b*c - 2*Re(conj(f)*d*e)
        
        DOUBLE p, sqrt_p, q, c, s, phi;
        
        p = SQR(m) - DOUBLE(3.0)*c1;
        q = m*(p - (DOUBLE(3.0)/DOUBLE(2.0))*c1) - (DOUBLE(27.0)/DOUBLE(2.0))*c0;
        sqrt_p = sqrt(DABS(p));
        
        phi = DOUBLE(27.0) * ( DOUBLE(0.25)*SQR(c1)*(p - c1) + c0*(q + DOUBLE(27.0)/DOUBLE(4.0)*c0));
        phi = (DOUBLE(1.0)/DOUBLE(3.0)) * atan2(sqrt(DABS(phi)), q);
        
        c = sqrt_p*cos(phi);
        s = (DOUBLE(1.0)/D_SQRT3)*sqrt_p*sin(phi);
        
        w[1]  = (DOUBLE(1.0)/DOUBLE(3.0))*(m - c);
        w[2]  = w[1] + s;
        w[0]  = w[1] + c;
        w[1] -= s;
        
        return 1;
    }
    
    // ----------------------------------------------------------------------------
    // Reduces a hermitian 3x3 matrix to real tridiagonal form by applying
    // (unitary) Householder transformations:
    //            [ d[0]  e[0]       ]
    //    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
    //            [       e[1]  d[2] ]
    // The function accesses only the diagonal and upper triangular parts of
    // A. The access is read-only.
    // ---------------------------------------------------------------------------
    
    void HouseholderTransform(COMPLEX A[3][3],COMPLEX Q[3][3],DOUBLE d[3], DOUBLE e[2]){
        
        COMPLEX u[3],q[3];
        COMPLEX omega,f;
        DOUBLE K, h, g;
        
        // Initialize Q to the identitity matrix
        for (int i=0; i < 3; i++){
            Q[i][i] = DOUBLE(1.0);
            for (int j=0; j < i; j++){
                Q[i][j] = Q[j][i] = DOUBLE(0.0);
            }
        }
        
        // Bring first row and column to the desired form
        h = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
        
        if (real(A[0][1]) > 0){
            g = -sqrt(h);
        }
        
        else{
            g = sqrt(h);
        }
        
        e[0] = g;
        f    = g * A[0][1];
        u[1] = conj(A[0][1]) - g;
        u[2] = conj(A[0][2]);
        
        omega = h - f;
        if (real(omega) > DOUBLE(0.0)){
            
            omega = DOUBLE(0.5) * (DOUBLE(1.0) + conj(omega)/omega) / real(omega);
            K = DOUBLE(0.0);
            for (int i=1; i < 3; i++){
                f    = conj(A[1][i]) * u[1] + A[i][2] * u[2];
                q[i] = omega * f;                  // p
                K   += real(conj(u[i]) * f);      // u* A u
            }
            K *= DOUBLE(0.5) * SQR_ABS(omega);
            
            for (int i=1; i < 3; i++){
                q[i] = q[i] - K * u[i];
            }
            
            d[0] = real(A[0][0]);
            d[1] = real(A[1][1]) - DOUBLE(2.0)*real(q[1]*conj(u[1]));
            d[2] = real(A[2][2]) - DOUBLE(2.0)*real(q[2]*conj(u[2]));
            
            // Store inverse Householder transformation in Q
            for (int j=1; j < 3; j++){
                f = omega * conj(u[j]);
                for (int i=1; i < 3; i++){
                    Q[i][j] = Q[i][j] - f*u[i];
                }
            }
            
            // Calculate updated A[1][2] and store it in f
            f = A[1][2] - q[1]*conj(u[2]) - u[1]*conj(q[2]);
        }
        
        else{
            
            for (int i=0; i < 3; i++){
                d[i] = real(A[i][i]);
            }
            
            f = A[1][2];
        }
        
        // Make (23) element real
        e[1] = abs(f);
        
        if (e[1] != DOUBLE(0.0)){
            
            f = conj(f) / e[1];
            
            for (int i=1; i < 3; i++){
                Q[i][2] = Q[i][2] * f;
            }
        }
    }
    
    // ----------------------------------------------------------------------------
    // Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
    // matrix A using the QL algorithm with implicit shifts, preceded by a
    // Householder reduction to real tridiagonal form.
    // The function accesses only the diagonal and upper triangular parts of A.
    // The access is read-only.
    // ----------------------------------------------------------------------------
    // Parameters:
    //   A: The hermitian input matrix
    //   Q: Storage buffer for eigenvectors
    //   w: Storage buffer for eigenvalues
    // ----------------------------------------------------------------------------
    // Return value:
    //   0: Success
    //  -1: Error (no convergence)
    // ----------------------------------------------------------------------------
    
    int EigensystemQL(COMPLEX A[3][3], COMPLEX Q[3][3], DOUBLE w[3]){
        
        DOUBLE e[3];                 // The third element is used only as temporary workspace
        DOUBLE g, r, p, f, b, s, c;  // Intermediate storage
        COMPLEX t;
        int m;
        
        //ITERATION COUNTER
        int nIter;
        int MaxIter=50;
        
        // Transform A to real tridiagonal form by the Householder method
        HouseholderTransform(A, Q, w, e);
        
        // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
        // with the QL method
        
        // Loop over all off-diagonal elements
        for (int l=0; l < 2; l++){
            
            nIter = 0;
            while (1){
                
                // Check for convergence and exit iteration loop if off-diagonal
                // element e(l) is zero
                for (m=l; m <= 1; m++){
                    
                    g = DABS(w[m])+DABS(w[m+1]);
                    
                    if (DABS(e[m]) + g == g){
                        break;
                    }
                }
                
                if (m == l){break;}
                
                //STOP AFTER 50 ITERATIONS IF NO CONVERGENCE CAN BE ACHIEVED
                if(nIter++ >= MaxIter){return -1;}
                
                // Calculate g = d_m - k
                g = (w[l+1] - w[l]) / (e[l] + e[l]);
                r = sqrt(SQR(g) + DOUBLE(1.0));
                
                if (g > 0){
                    g = w[m] - w[l] + e[l]/(g + r);
                }
                
                else{
                    g = w[m] - w[l] + e[l]/(g - r);
                }
                
                s = c = DOUBLE(1.0);
                p = DOUBLE(0.0);
                
                for (int i=m-1; i >= l; i--){
                    
                    f = s * e[i];
                    b = c * e[i];
                    
                    if (DABS(f) > DABS(g)){
                        
                        c      = g / f;
                        r      = sqrt(SQR(c) + DOUBLE(1.0));
                        e[i+1] = f * r;
                        c     *= (s = DOUBLE(1.0)/r);
                    }
                    
                    else{
                        
                        s      = f / g;
                        r      = sqrt(SQR(s) + DOUBLE(1.0));
                        e[i+1] = g * r;
                        s     *= (c = DOUBLE(1.0)/r);
                        
                    }
                    
                    g = w[i+1] - p;
                    r = (w[i] - g)*s + DOUBLE(2.0)*c*b;
                    p = s * r;
                    w[i+1] = g + p;
                    g = c*r - b;
                    
                    // Form eigenvectors
                    for (int k=0; k < 3; k++){
                        t = Q[k][i+1];
                        Q[k][i+1] = s*Q[k][i] + c*t;
                        Q[k][i]   = c*Q[k][i] - s*t;
                    }
                    
                }
                
                w[l] -= p;
                e[l]  = g;
                e[m]  = DOUBLE(0.0);
            }
        }
        
        return 0;
    }
    
    
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // Calculates the eigenvalues and Normalized eigenvectors of a hermitian 3x3
    // matrix A using Cardano's method for the eigenvalues and an analytical
    // method based on vector cross products for the eigenvectors. However,
    // if conditions are such that a large Error in the results is to be
    // expected, the routine falls back to using the slower, but more
    // accurate QL algorithm. Only the diagonal and upper triangular parts of A need
    // to contain meaningful values. Access to A is read-only.
    // ----------------------------------------------------------------------------
    // Parameters:
    //   A: The hermitian input matrix
    //   Q: Storage buffer for eigenvectors
    //   w: Storage buffer for eigenvalues
    // ----------------------------------------------------------------------------
    // Return value:
    //   0: Success
    //  -1: Error
    // ----------------------------------------------------------------------------
    
    int Eigensystem(COMPLEX A[3][3], COMPLEX Q[3][3], DOUBLE w[3]){
        
        DOUBLE Norm;          // Squared Norm or inverse Norm of current eigenvector
        DOUBLE Error;         // Estimated maximum roundoff Error
        DOUBLE MaxAbsEval,Estimator;          // Intermediate storage
        
        // Calculate eigenvalues
        Eigenvalues(A, w);
        
        //Calculate estimated Error
        MaxAbsEval=std::max(DABS(w[0]),DABS(w[1]));
        MaxAbsEval=std::max(MaxAbsEval,DABS(w[2]));
        
        if(MaxAbsEval<DOUBLE(1.0)){
            Estimator=MaxAbsEval;
        }
        else{
            Estimator=SQR(MaxAbsEval);
        }
        
        Error = DOUBLE(256.0) * DBL_EPSILON * SQR(Estimator);
        
        Q[0][1] = A[0][1]*A[1][2] - A[0][2]*real(A[1][1]);
        Q[1][1] = A[0][2]*conj(A[0][1]) - A[1][2]*real(A[0][0]);
        Q[2][1] = SQR_ABS(A[0][1]);
        
        // Calculate first eigenvector by the formula
        //   v[0] = conj( (A - w[0]).e1 x (A - w[0]).e2 )
        Q[0][0] = Q[0][1] + A[0][2]*w[0];
        Q[1][0] = Q[1][1] + A[1][2]*w[0];
        Q[2][0] = (real(A[0][0]) - w[0]) * (real(A[1][1]) - w[0]) - Q[2][1];
        
        Norm    = SQR_ABS(Q[0][0]) + SQR_ABS(Q[1][0]) + SQR(real(Q[2][0]));
        
        // If vectors are nearly linearly dependent, or if there might have
        // been large cancellations in the calculation of A(I,I) - W(1), fall
        // back to QL algorithm.
        
        // Note that this simultaneously ensures that multiple eigenvalues do
        // not cause problems: If W(1) = W(2), then A - W(1) * I has rank 1,
        // i.e. all columns of A - W(1) * I are linearly dependent.
        
        //CHECK ACCEPTANCE CRITERION AND NORMALIZE
        if(Norm>Error){
            Norm = sqrt(DOUBLE(1.0) / Norm);    
            for (int j=0; j < 3; j++){
                Q[j][0] = Q[j][0] * Norm;
            }
        }
        
        //SWITCH TO QR ALGORITHM IF REJECTED
        else{
            return EigensystemQL(A, Q, w);
        }
        
        // Calculate second eigenvector by the formula
        //   v[1] = conj( (A - w[1]).e1 x (A - w[1]).e2 )
        Q[0][1]  = Q[0][1] + A[0][2]*w[1];
        Q[1][1]  = Q[1][1] + A[1][2]*w[1];
        Q[2][1]  = (real(A[0][0]) - w[1]) * (real(A[1][1]) - w[1]) - real(Q[2][1]);
        
        Norm     = SQR_ABS(Q[0][1]) + SQR_ABS(Q[1][1]) + SQR(real(Q[2][1]));
        
        //CHECK ACCEPTANCE CRITERION AND NORMALIZE
        if(Norm>Error){
            
            Norm = sqrt(DOUBLE(1.0) / Norm);
            
            for (int j=0; j < 3; j++){
                Q[j][1] = Q[j][1] * Norm;
            }
        }
        
        //SWITCH TO QR ALGORITHM IF REJECTED
        else{
            return EigensystemQL(A, Q, w);
        }
        
        // Calculate third eigenvector according to
        //   v[2] = conj(v[0] x v[1])
        Q[0][2] = conj(Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1]);
        Q[1][2] = conj(Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1]);
        Q[2][2] = conj(Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1]);
        
        return 1;
    }
    
}

#endif
