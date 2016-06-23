#ifndef __DEFINITIONS_CPP__
#define __DEFINITIONS_CPP__

#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

//DEFINITIONS OF DATATYPES
typedef int INT;
typedef double DOUBLE;
typedef std::complex<DOUBLE> COMPLEX;

//DETERMINE MAXIMUM WORKING PRECISION
static const int MAX_DIGITS_PRECISION=std::numeric_limits<DOUBLE>::digits10;
static const int OUTPUT_PRECISION=MAX_DIGITS_PRECISION;

using std::sqrt;
using std::pow;
using std::abs;
using std::cos;
using std::sin;
using std::exp;
using std::atan2;

//BASIC CONSTANTS
#ifndef ComplexI
#define ComplexI COMPLEX(0.0,1.0) // i
#endif

#ifndef PI
#define PI         DOUBLE(3.141592653589793238462643383279502884197169399375105820974944592)   // pi
#endif

#define D_SQRT2    DOUBLE(1.414213562373095048801688724209698078569671875376948073176679738)   // sqrt(2)
#define D_SQRT3    DOUBLE(1.732050807568877293527446341505872366942805253810380628055806979)   // sqrt(3)

//BASIC MACROS
#define DABS(x)  (((x)>=(0))?(x):(-x))  // |x|
#define SQR(x)      ((x)*(x))                        // x^2 
#define SQR_ABS(x)  (SQR(real(x)) + SQR(imag(x)))  // |x|^2
#define MOD(x,y)    ((x)>=(0)?((x)%(y)):(((x)+(y))%(y)))       // x mod y
#define DELTA(x,y)    (((x)==(y))?(1):(0))       // Kronecker delta
#define DSQRT(x) (((x)>=0)?(sqrt(x)):(ComplexI*sqrt(-(x)))) // PRINCIPAL SQUARE ROOT OF REAL NUMBER
#define SIGN(x)  ((x)>0?(1):((x)<0?(-1):(0))) // SIGN FUNCTION
#define LEVICIVITA4D(i,j,k,l) (((i)-(j))*((i)-(k))*((i)-(l))*((j)-(k))*((j)-(l))*((k)-(l))/DOUBLE(12.0)) // 4D e_{ijkl}


//COMPILER FLAGS FOR STATIC BOX OR EXPANDING BOX
#define MINKOWSKI_FLAG 100
#define BJORKEN_FLAG 101

//COMPILER FLAGS FOR SUNc GAUGE GROUP
#define U1_FLAG  111
#define SU2_FLAG 222
#define SU3_FLAG 333

//COMPILER FLAGS INITIAL CONDITIONS
#define THERMAL_FLAG 19
#define QP_FLAG 91
#define LOAD_FLAG 666

// MPI //
#include <mpi.h>

// OPEN MP //
#include <omp.h>

// FFTW INCLUSION //
#include "fftw3.h"

#ifndef MY_FFTW_PLANNER_FLAG
#define MY_FFTW_PLANNER_FLAG FFTW_MEASURE
#endif

#endif
