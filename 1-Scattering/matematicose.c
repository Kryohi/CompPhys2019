

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

inline double variance(const double * A, size_t length)
{
    double * A2 = malloc(length * sizeof(double));
    elforel(A,A,A2,length);
    double var = mean(A2,length) - mean(A,length)*mean(A,length);
    free(A2);
    return var;
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


 
