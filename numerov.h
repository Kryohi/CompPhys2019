
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <complex.h>
//#include <fftw3.h>
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


Spectrum numerov(int nmax, int l, int xmax, double rmax, double h2m, double Estep, bool normalize, dArray bc0, double (*f)(double));
void numerov_forward(double h, int xc, int xmin, const double * k2, double * y);
void numerov_backward(double h, int xc, int xmax, const double * k2, double * y);
double V_ho(double x);
double E0(double (*V)(double), double h, double rmax);
double E0_stupid(double (*V)(double), double h, double rmax);
double normalizationFactor(const double * eigv, double h, int x1, int x2);
int nodeNumber(const double * eigv, size_t N);
void save2csv(Spectrum * spectra, int lmax, int nmax, int xmax);

