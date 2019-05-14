#include "hydrogen.h"

#define eigen_method_diy 0
#define h2m 0.5

/* STO-3G parameters */
#define a1 0.109818
#define a2 0.405771
#define a3 2.22726


int main(int argc, char** argv)
{ 
    
    multiplytest();
    double a[] = {1,0,0,1};
    double b[] = {2,0,0,2};
    double *eval = malloc(2 * sizeof(double));
    eval = gen_eigenvalues_nxn(a, b, 2);
    printf ("[ %g, %g\n", eval[0], eval[1]);
    printf ("  %g, %g ]\n", eval[2], eval[3]);
    
    return 0;
}

double * gen_eigenvalues_nxn(double *a, double *b, uint16_t n)
{
    double *eigenvalues = malloc(n * sizeof(double));
    double *eigenvectors = malloc(n*n * sizeof(double));
    zeros(n, eigenvalues); zeros(n*n, eigenvectors);
    
    gsl_matrix_view A = gsl_matrix_view_array(a, n, n);
    gsl_matrix_view B = gsl_matrix_view_array(b, n, n);
    gsl_vector_view eval = gsl_vector_view_array(eigenvalues, n);
    gsl_matrix_view evec = gsl_matrix_view_array(eigenvectors, n, n);
    
    
    #if eigen_method_diy
    /* Metodo 1:
     * algoritmo di Boh
    */
    
    
    #else
    /* Metodo 2:
     * gsl_eigen_gensymmv()
    */
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(&A.matrix, &B.matrix, &eval.vector, &evec.matrix, w);
    gsl_eigen_gensymmv_free(w);

    #endif
    
    
    return eigenvalues;
}

void multiplytest(void) {
    double a[] = { 0.11, 0.12, 0.13,
                 0.21, 0.22, 0.23 };

    double b[] = { 1011, 1012,
                 1021, 1022,
                 1031, 1032 };

    double c[] = { 0.00, 0.00,
                 0.00, 0.00 };

    gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
    gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
    gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);

    /* Compute C = A B */

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

    printf ("[ %g, %g\n", c[0], c[1]);
    printf ("  %g, %g ]\n", c[2], c[3]);

}
