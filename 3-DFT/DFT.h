#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <complex.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_permutation.h>
#include <fftw3.h>
#include "../matematicose.c"
#include "../misccose.c"
#include "../numerov.c"

double local_energy(double rho);

// Kohn-Sham potential terms
double V_ks(double r, const double *rho);
double V_ks_rho(double r); // same as above, but uses the global variable RHO
double V_ext(double r);
double V_h(double r, const double *rho);

// Kohn-Sham energy functional terms
double E_ks(double *rho);
double E_ext(double *rho);
double E_H(double *rho);
double E_XC(double *rho);
double T_S(double *rho);

