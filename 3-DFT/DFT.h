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

double V_ks(double r, double *rho);
double local_energy(double rho);
double E_ks(double *rho);


