#pragma once
#include <math.h>
//#include <stdatomic.h>
//#include <pthread.h>


void pointwise(double (*f)(double), double * A, size_t length);
int double_max_index(double * A, size_t length);
int double_min_index(double * A, size_t length);
double sum(const double *A, size_t length);
int intsum(const int * A, size_t length);
double mean(const double * A, size_t length);
double intmean(const int * A, size_t length);
double variance(const double * A, size_t length);
double variance_corr(const double * A, double tau, size_t length);
void zeros(size_t length, double *A);
void elforel(const double *A, const double * B, double * C, size_t length);
bool isApproxEqual(double a, double b);
double zerosecant(double (*f)(double), double x1, double x2, double inf, double sup);
double secant(double (*f)(double), double c, double x1, double x2, double inf, double sup);
double findzero_last(double (*f)(double), double c, double x1, double x2, double inf, double sup);
void fast_bessel(double x, double lmax, double * J);
double der3(double * F, int x, double h);
double der5(double * F, int x, double h);
double der5_c(double (*f)(double), double x, double h);
double simpson_integral(double *fun, int xmax, double h);
void simple_acf(const double *H, size_t length, int k_max, double * acf);
void fft_acf(const double *H, size_t length, int k_max, double * acf);
double grad_descent_1D(double (*f)(double), double x1, double x2);
double stochastic_grad_descent_1D(double (*f)(double), double x1, double x2);