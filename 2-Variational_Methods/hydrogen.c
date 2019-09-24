#include "hydrogen.h"

#define eigen_method_diy 0
#define h2m 0.5

/* STO-3G parameters */
#define a1 0.109818
#define a2 0.405771
#define a3 2.22726


void multiplytest(void);


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

void multiplytest(){
    double a[] = { 0.11, 0.12, 0.13,
                 0.21, 0.22, 0.23,
                 0.31, 0.32, 0.33};

    double b[] = { 1011, 1012, 1013,
                 1021, 1022, 1023,
                 1031, 1032, 1033 };

    double c[] = { 1.00, 0.00, 0.00,
                 0.00, 2.00, 0.00,
                 0.00, 0.00, 3.00};
                 
    double d[] = { 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00};

    double e[] = { 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00};
                 
    gsl_matrix_view A = gsl_matrix_view_array(a, 3, 3);
    gsl_matrix_view B = gsl_matrix_view_array(b, 3, 3);
    gsl_matrix_view C = gsl_matrix_view_array(c, 3, 3);
    gsl_matrix_view D = gsl_matrix_view_array(d, 3, 3);
    gsl_matrix_view E = gsl_matrix_view_array(e, 3, 3);
    // we want to compute ACB, where C is the diagonal matrix. To do this we fist compute D=CB

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &C.matrix, &B.matrix,
                  0.0, &D.matrix);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &D.matrix,
                  0.0, &E.matrix);
    
    printf ("[ %g, %g, %g \n", e[0], e[1], e[2]);
    printf ("  %g, %g, %g \n", e[3], e[4], e[5]);
    printf ("  %g, %g, %g]\n", e[6], e[7], e[8]);
}
