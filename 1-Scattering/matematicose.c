

/*
    Simple math
*/

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
double zerosecant(double (*f)(double), double x1, double x2, double inf, double sup);
double secant(double (*f)(double), double c, double x1, double x2, double inf, double sup);
double findzero_last(double (*f)(double), double c, double x1, double x2, double inf, double sup);
void fast_bessel(double x, double lmax, double * J);
double der3(double * F, int x, double h);
double der5(double * F, int x, double h);
double der5_c(double (*f)(double), double x, double h);
double integrale_simpson(double *fun, int xmax, double h);

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
    double step = (x2-x1)/20;
    for (double x=x2; x>x1; x -= step)
        if ((f(x)-c)*(f(x-step)-c) < 0)
            return secant(f, c, x-step, x, inf, sup);
    perror("no zeros found");
    return -1;
}

void fast_bessel(double x, double lmax, double * J)
{
    //double j_1 = cos(x)/x;
    J[0] = sin(x)/x;
    J[1] = (J[0] - cos(x))/x;
    
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

// numerical derivatives

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


double integrale_simpson(double *fun, int xmax, double h){
    //xmax should be even, given the algorithm used, and the value of the vector should be calculated equally spaced by h
    double integral = 0.;

    for(int i = 1; i < xmax-1; i+=2) {
        integral += h*(fun[i-1] + 4.*fun[i]+fun[i+1])/3.;
    }
    
    return integral;
}


