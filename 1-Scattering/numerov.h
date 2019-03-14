 

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
    size_t length;
    double *data;
} Spectrum;


void numerov_forward(int l, double E, double h, double rmax, double * x);
double E0(double h, double rmax);


