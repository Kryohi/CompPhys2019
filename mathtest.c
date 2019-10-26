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
#include "matematicose.c"
#include "misccose.c"

double test_function1(double * x);

int main(int argc, char** argv)
{ 
    int d = 2;
    double h = 1e-4;
    double y[d];
    double x[] = {M_PI/6, M_PI/6};
    
    // Test funzione gradiente
    for (int i=0; i<d; i++) {
        y[i] = der5_part(test_function1, x, d, i, h);
        printf("f(x) = %f\n", test_function1(x));
        printf("df/d%d = %f\n", i, y[i]);
    }
    
    // Test gradient descent
    double minimum[d];
    grad_descent(test_function1, -2.0, 3.0, d, minimum);
    
    printf("\nxmin = %f, %f", minimum[0], minimum[1]);
    
    return 0;
}

double test_function1(double * x)
{
    return sin(x[0]) + cos(x[1]);
}



