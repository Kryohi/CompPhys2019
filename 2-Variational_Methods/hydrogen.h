#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <fftw3.h>
#include "../matematicose.c"
#include "../misccose.c"

typedef struct Spectrum { 
    size_t xmax;
    int nmax;
    double * EE;
    double * eigfuns;  // xmax*nmax 1D array
} Spectrum;

typedef struct dArray { 
    size_t length;
    double * data;
} dArray;



double * eigenvalues_nxn(double *a, double *b, uint16_t n);


void multiply3matrix (double *a, double *c, double *b, double *e, uint16_t n){ //the fourth input is the result of the multiplication, 
    gsl_matrix_view A = gsl_matrix_view_array(a, n, n);
    gsl_matrix_view B = gsl_matrix_view_array(b, n, n);
    gsl_matrix_view C = gsl_matrix_view_array(c, n, n);
    
    double d[] = { 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00};


    gsl_matrix_view D = gsl_matrix_view_array(d, n, n);
    gsl_matrix_view E = gsl_matrix_view_array(e, n, n);
     // we want to compute ACB, where C is the diagonal matrix. To do this we fist compute D=CB
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &C.matrix, &B.matrix,
                  0.0, &D.matrix);
    // now we do E=AD 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &D.matrix,
                  0.0, &E.matrix);
    
}
