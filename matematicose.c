
#include "matematicose.h"

/*
    Simple math
*/


inline double sum(const double * A, size_t length)   
{
    double s = 0.;
    for (int i=0; i<length; i++)
        s += A[i];

    return s;
}

inline int intsum(const int * A, size_t length)   
{
    int s = 0;
    for (int i=0; i<length; i++)
        s += A[i];

    return s;
}

inline double dot(const double * A, double * B, size_t length)  
{
    double result = 0.0;

    for (int i=0; i<length; i++)
        result += A[i]*B[i];

    return result;
}

inline void elforel(const double * A, const double * B, double * C, size_t length)  
{
    for (int i=0; i<length; i++)
        C[i] = A[i]*B[i];
}

inline double mean(const double * A, size_t length) 
{
    return sum(A,length)/length;
}

inline double intmean(const int * A, size_t length) 
{
    int s = 0;
    for (int i=0; i<length; i++)
        s += A[i];
    
    return (double)s/length;
}

inline void zeros(size_t length, double *A)
{
    for (int i=length; i!=0; i--)
        A[i] = 0.0;
}

inline double min_d(const double * A, size_t length)
{
    double min = A[0];
    for (int i=1; i<length; i++)
        if (A[i] < min) min = A[i];
    
    return min;
}

inline double max_d(const double * A, size_t length)
{
    double max = A[0];
    for (int i=1; i<length; i++)
        if (A[i] > max) max = A[i];
    
    return max;
}

inline double min_absd(const double * A, size_t length)
{
    double min = fabs(A[0]);
    for (int i=1; i<length; i++)
        if (fabs(A[i]) < min) min = fabs(A[i]);
    
    return min;
}

inline double max_absd(const double * A, size_t length)
{
    double max = fabs(A[0]);
    for (int i=1; i<length; i++)
        if (fabs(A[i]) > max) max = fabs(A[i]);
    
    return max;
}


// apply a function to each element in the A array 
// TODO check the slowdown vs a single for 
void pointwise(double (*f)(double), double * A, size_t length) 
{
    for (int n=0; n<length; n++)
        A[n] = f(A[n]);
}
// apply a function with 2 parameters and store the results in an array B half the lenght of A
/*void pointwise2D(double (*f)(double *), double * A, double * B, size_t length) 
{
    for (int n=0; n<length; n++)
        B[n] = f(A[2*n], A[2*n+1]);
}*/


// return the index of the element with maximum value in a double array
int double_max_index(double * A, size_t length)
{
    int idx = 0;
    for (int n=1; n<length; n++)
        if (A[n]>A[idx]) idx = n;
    
    return idx;
}
int double_min_index(double * A, size_t length)
{
    int idx = 0;
    for (int n=1; n<length; n++)
        if (A[n]<A[idx]) idx = n;
    
    return idx;
}


inline double variance(const double * A, size_t length)
{
    double * A2 = malloc(length * sizeof(double));
    elforel(A,A,A2,length);
    double var = mean(A2,length) - mean(A,length)*mean(A,length);
    free(A2);
    return var;
}

// Secant method with a function as input
double zerosecant(double (*f)(double), double x1, double x2, double inf, double sup)
{
    if (f(x1)>inf && f(x1)<sup)
        return x1;
    else if (f(x2)>inf && f(x2)<sup)
        return x2;
    else if (f(x1)*f(x2)>0) {
        perror("f(X1) and f(X2) must have an opposing sign");
        return -1;
    }
    else
    {
        double x_;
        while (f(x2)<inf || f(x2)>sup)  // restituisce il valore dove la funzione fa circa 0
        {
            x_ = x2;  // x precedente
            x2 = x2 - f(x2)*(x2-x1)/(f(x2)-f(x1));
            x1 = x_;
        }
    }
    return x2;
}

// same as above, but takes an additional argument that subctracts from the function
// i.e. finds where f(x) and c meet
double secant(double (*f)(double), double c, double x1, double x2, double inf, double sup)
{
    //printf("V(x1) = %f\tV(x2) = %f\t c=%f\n", f(x1), f(x2), c);
    if ((f(x1)-c)>inf && (f(x1)-c)<sup)
        return x1;
    else if ((f(x2)-c)>inf && (f(x2)-c)<sup)
        return x2;
    else if ((f(x1)-c)*(f(x2)-c)>0) {
        perror("f(X1) and f(X2) must have an opposing sign");
        return -1;
    }
    else
    {
        double x_;
        while ((f(x2)-c)<inf || (f(x2)-c)>sup)  // restituisce il valore dove la funzione fa circa 0
        {
            x_ = x2;  // x precedente
            x2 = x2 - (f(x2)-c)*(x2-x1)/((f(x2)-c)-(f(x1)-c));
            x1 = x_;
            //printf("x2 fin = %f\n",x2);
        }
        
    }
    return x2;
}

double findzero_last(double (*f)(double), double c, double x1, double x2, double inf, double sup)
{
    double step = (x2-x1)/1000;
    for (double x=x2; x>x1; x -= step)  {
        //printf("%f\t%f\n", f(x)-c, f(x-step)-c);
        if ((f(x)-c)*(f(x-step)-c) < 0) {
            //printf("\ninversion at %f - %f\n", x-step, x);
            return secant(f, c, x-step, x, inf, sup);}
    }
    perror("no zeros found");
    return -1;
}

// J should already provide the 3 starting points
inline void fast_bessel(double x, double lmax, double * J)
{    
    for (int l=1; l<lmax; l++)
        J[l+1] = ((2*l+1)/x) * J[l] - J[l-1];
}


/*
 * Put in the array A gaussian-distributed numbers around 0, with standard deviation sigma
*/

inline void vecBoxMuller(double sigma, size_t length, double * A)
{
    double x1, x2;

    for (int i=0; i<round(length/2); i++) {
        x1 = (double) rand() / (RAND_MAX + 1.0);
        x2 = (double) rand() / (RAND_MAX + 1.0);
        A[2*i] = sigma * sqrt(-2*log(1-x1)) * cos(2*M_PI*x2);
        A[2*i+1] = sigma * sqrt(-2*log(1-x2)) * sin(2*M_PI*x1);
    }
}



/*
 * 
 * Calculus
 * 
*/

// Numerical derivatives

double der3(double * F, int x, double h)
{
    return (F[x+1] - F[x-1])/(2*h);
}

double der5(double * F, int x, double h)
{
    return (-F[x+2] + 8*F[x+1] - 8*F[x-1] + F[x-2])/(12*h);
}

double der5_c(double (*f)(double), double x, double h)
{
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(12*h);
}

double der5_part(double (*f)(double*), double * x, size_t d, int l, double h)
{
    double xh[d], x2h[d], x_h[d], x_2h[d];
    
    for (int i=0; i<d; i++)
    {
        xh[i] = x[i];
        x2h[i] = x[i];
        x_h[i] = x[i];
        x_2h[i] = x[i];
    }
    xh[l] = x[l]+h;
    x2h[l] = x[l]+2*h;
    x_h[l] = x[l]-h;
    x_2h[l] = x[l]-2*h;
    
    return (-f(x2h) + 8*f(xh) - 8*f(x_h) + f(x_2h))/(12*h);
}


// Numerical integration

double simpson_integral(double *fun, int xmax, double h){
    //xmax should be even, given the algorithm used, and the value of the vector should be calculated equally spaced by h
    double integral = 0.;

    for(int i = 1; i < xmax-1; i+=2) {
        integral += h*(fun[i-1] + 4.*fun[i]+fun[i+1])/3.;
    }
    
    return integral;
}




// Gradient descent (toward madness)

double grad_descent_1D(double (*f)(double), double x1, double x2)
{
    double x = (x2-x1)/2; // starting point
    double scale = fabs(f(x2)-f(x)); // to get an order of magnitude of the variation of V
    double gamma = scale/200; // learning rate
    double h = (x2-x1)/5e4;
    double grad = 10.;
    
    while (fabs(grad) > 1e-7)
    {
        grad = der5_c(f,x,h);
        x -= gamma*grad;
    }
    //printf("grad = %0.10f\tr=%f\n", grad, r);
    
    return x;
}



// First iterates through some randomly placed points X0, in the region bounded by x1 and x2,
// in order to avoid regions too far from the minimum that might have vanishing gradient.
// Then starts the gradient descent from the point with the lowest value found.
// The resulting minimum (array) is put into x

void grad_descent(double (*f)(double*), double x1, double x2, uint16_t d, double * x, bool ascent)
{
    srand(42);
    int N = 50*d; // number of different starting points, the linear dependence on d was chosen arbitrarily
    
    // create grid of N randomly placed points in a hypercube [-x1, x2]^d
    //TODO: use Poisson-disk sampling or another pseudorandom algorithm instead of simple rand()
    double X0[N][d];
    double F0[N];
    for (int n=0; n<N; n++) {
        for (int i=0; i<d; i++) {
            X0[n][i] = (rand()/(RAND_MAX+1.))*fabs(x2-x1) + x1;
        }
    }
    
    double temp_X[d];
    for (int n=0; n<N; n++) {
        for (int i=0; i<d; i++)
            temp_X[i] = X0[n][i];
            
        //F0[n] = f(X0[n]);
        F0[n] = f(temp_X);
    }
    
    int idxmin = double_min_index(F0, N);
    int idxmax = double_max_index(F0, N);
    if (ascent) {
        for (int i=0; i<d; i++)
            x[i] = X0[idxmax][i];   // starting point for the gradient ascent
    }
    else    {
        for (int i=0; i<d; i++)
            x[i] = X0[idxmin][i];   // starting point for the gradient descent
    }
    
    // to get an order of magnitude of the variation of V
    double scale = fabs(F0[double_max_index(F0, N)] - F0[double_min_index(F0, N)]);
    printf("\nscale = %f - %f = %f", F0[double_max_index(F0, N)], F0[double_min_index(F0, N)], scale);
    double gamma = 0.00001; // learning rate
    double h = fabs(x2-x1)/1e6;
    double grad[d];
    for (int i=0; i<d; i++)
        grad[i] = 10.0;
    
    int count = 0;
    int sign;
    if (ascent) {
        sign = -1;
    }
    else    {
        sign = 1;
    }
    
    while (max_absd(grad,d) > 1e-4)
    {
        for (int i=0; i<d; i++)
        {
            grad[i] = der5_part(f,x,d,i,h); // spostare qui invece di usare funzione?
            x[i] -= sign*gamma*grad[i];
            printf("df/da%d = %f\n", i, grad[i]);
            printf("a[%d] = %f\n", i, x[i]);
            printf("\n");
        }
        count++;
    }
    printf("\nxmin = %f, %f\t %d steps performed", x[0], x[1], count);
    
}



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



