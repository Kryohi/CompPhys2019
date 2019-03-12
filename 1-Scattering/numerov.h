 

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


void numerov(int n, int l, double h, double rmax, double * x)

void simple_acf(const double *H, size_t length, int k_max, double * acf);
void fft_acf(const double *H, size_t length, int k_max, double * acf);
double sum(const double *A, size_t length);
int intsum(const int * A, size_t length);
double mean(const double * A, size_t length);
double intmean(const int * A, size_t length);
double variance(const double * A, size_t length);
double variance_corr(const double * A, double tau, size_t length);
void zeros(size_t length, double *A);
void elforel(const double *A, const double * B, double * C, size_t length);
bool isApproxEqual(double a, double b);


