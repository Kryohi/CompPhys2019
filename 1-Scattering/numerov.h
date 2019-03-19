
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <complex.h>
#include <fftw3.h>


typedef struct DoubleArray { 
    size_t length;
    double *data;
} DoubleArray;

typedef struct Spectrum { 
    size_t xmax;
    int nmax;
    double * EE;
    double * eigfuns;  // xmax*nmax 1D array
} Spectrum;


void numerov_forward(double h, int xc, const double * k2, double * y);
void numerov_backward(double h, int xc, int xmax, const double * k2, double * y);
Spectrum numerov(int nmax, int l, int xmax, double rmax, double Estep, double (*f)(double));
double V_ho(double x);
double E0(double h, double rmax, double (*V)(double));


