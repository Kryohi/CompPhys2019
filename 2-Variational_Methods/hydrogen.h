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


