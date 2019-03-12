
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>


#define HBAR_SPERIMENTALI 1.05457180013e-39

// number 
#define LMAXX 10

// number 
#define NMAXX 10

// number 
#define MAX_MESH_POINTS 200


// Harmonic potential
inline double V(double r)
{
    return 0.5*r*r;
}

// returns the (partial) array of eigenvalues and a (partial) matrix of eigenvectors
void Hsolve(int n, int L, int N, double hbarm; double rmax, double E[N], double u[N][n])
{
    int i,j;
    double d, e1, e2, h, r, vv;
    
    gsl_vector *eval = gsl_vector_calloc(n);     // eigenvalues array
    gsl_matrix *evec = gsl_matrix_calloc(n,n);   // eigenvectors array
    gsl_matrix *A    = gsl_matrix_calloc(n,n);   // my big fat matrix
    
    // step length
    h = rmax/n;
    
    // build up the tridiagonal matrix
    d = hbarm*2./(h*h);
    e1 = -hbarm/(h*h);
    
    for(i=0; i<n; i++)
    {
        r = (i+1)*h;
        vv = V(r);
        if (i<n-1) gsl_matrix_set(A,i,i, d+vv+hbarm*L*(L+1)/(r*r));
        if (i<n-2) gsl_matrix_set(A,i,i+1,e1);
        if (i>0)   gsl_matrix_set(A,i,i-1,e1);
    }
    // print matrix to terminal
    gsl_matrix_fprintf(stdout, A, "%f");
    
    // allocate temp memory to work with
    gsl_matrix_symmv_workspace *temp = gsl_eigen_symmv_alloc(n);
    
    // Diagonalization
    gsl_eigen_symmv(A, eval, evec, w);
    
    // Sorting of the eigenvalues
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    
    // copy the results in the returned array and matrix
    for(i=0; i<N; i++)
    {
        E[i] = gsl_vector_get(eval,i);
        for (j=0; j<n; j++)
            u[i][j] = gsl_matrix_get(evec,j,i);
    }
    
    //releases memory
    gsl_matrix_symmv_free(temp);
    gsl_matrix_free(A);
}


int main(int argc, char** argv)
{
    FILE *out;
    out = fopen("u.out", "w+");

    printf("Enter ");
    
    
    for (l=1; l<LMAX; l++)
    {
        Hsolve();
        
    }
    
    
    return 0;
}




