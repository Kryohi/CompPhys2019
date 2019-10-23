#include "hydrogen.h"


#define h2m 0.5

/* STO-3G parameters */
#define a1 0.109818
#define a2 0.405771
#define a3 0.222776 //NB PROBLEMA QUI: PEDERIVA SCRIVE 2.22776 ma altri ci mettono lo 0 davanti

// STO-2G parameters
#define a12g 0.151623
#define a22g 0.384244 //funziona meglio con gli altri ma non mi è chiaro perchè


void multiplytest(void);


// returns first eigenvalue
double gen_eigenvalues_nxn_gsl(double *a, double *b, uint16_t n)
{
    
    gsl_matrix_view A = gsl_matrix_view_array(a, n, n);
    gsl_matrix_view B = gsl_matrix_view_array(b, n, n);

    gsl_vector *eval1 = gsl_vector_alloc (n);
    gsl_matrix *evec1 = gsl_matrix_alloc(n, n);
    
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(&A.matrix, &B.matrix, eval1, evec1, w);
    gsl_eigen_gensymmv_free(w);

    for (int i = 0; i < n; i++)
      {
        double eval_i = gsl_vector_get (eval1, i);
        gsl_vector_view evec_i = gsl_matrix_column (evec1, i);

        //printf ("eigenvector = \n");
        //gsl_vector_fprintf (stdout,&evec_i.vector, "%g");
      }
      
    return gsl_vector_get(eval1, 0);
}


void gen_eigenvalues_nxn_man(double *a, double *b, uint16_t n){

    int i,j;

    gsl_matrix_view A = gsl_matrix_view_array(a, n, n);
    gsl_matrix_view B = gsl_matrix_view_array(b, n, n);
    
    gsl_vector *evalb = gsl_vector_alloc (n);
    gsl_matrix *evecb = gsl_matrix_alloc(n, n);
    
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(&B.matrix, evalb, evecb, w);
    gsl_eigen_symmv_free(w);
    
    double lambda1=sqrt(1./(gsl_vector_get (evalb, 0)));
    double lambda2=sqrt(1./(gsl_vector_get (evalb, 1)));
    double lambda3=sqrt(1./(gsl_vector_get (evalb, 2)));
    
    double lambdab[]={lambda1,0.,0.,  //matrix of eigenvalue of B inverted and square rooted
                     0.,lambda2,0.,
                     0.,0.,lambda3};
    
     double Phibp[] = { 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00};
                                 
    gsl_matrix_view lambdabm = gsl_matrix_view_array(lambdab, n, n);
    gsl_matrix_view Phibm = gsl_matrix_view_array(Phibp, n, n);  
    //gsl_matrix_view Evecb = gsl_matrix_view_array(evecb, n, n);   
                     
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, evecb, &lambdabm.matrix,
                  0.0, &Phibm.matrix);                 
    

    
    double Ap[] = { 0.00, 0.00, 0.00,
                    0.00, 0.00, 0.00,
                    0.00, 0.00, 0.00};
    
    double Phibt[] = { 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00,
                 0.00, 0.00, 0.00};

    gsl_matrix_view Phibtm = gsl_matrix_view_array(Phibt, n, n);

    gsl_matrix_transpose_memcpy(&Phibtm.matrix, &Phibm.matrix);
    
    multiply3matrix(Phibt, a, Phibp, Ap,n);
    
    gsl_matrix_view Apm= gsl_matrix_view_array(Ap, n, n);
    
    gsl_vector *evalAp = gsl_vector_alloc (n);
    gsl_matrix *evecAp = gsl_matrix_alloc(n, n);
    
    gsl_eigen_symmv_workspace * ww = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(&Apm.matrix, evalAp, evecAp, ww);
    gsl_eigen_symmv_free(ww);
    
    double Phi[]={ 0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00};
                   
    double evecbv[]={ 0.00, 0.00, 0.00, //stupida cosa per ottenere degli array da dare in pasto alla fn moltiplica matrici
                   0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00};
     int f=0;              
     for (i = 0; i < n; i++)
      {
        for (j=0; j<n;j++){
        
            evecbv[j+f]=gsl_matrix_get(evecb, i, j);
            f+=3;
        }
      }              
     double evecApv[]={ 0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00};
     int ff=0;              
     for (i = 0; i < n; i++)
      {
        for (j=0; j<n;j++){
            evecApv[j+ff]=gsl_matrix_get(evecAp, i, j);
            //printf("eigenstate of A=%g\n", evecApv[j+ff]);
            ff+=3;
        }
      }  
                   
    multiply3matrix(evecbv,lambdab,evecApv, Phi, n);
    
    gsl_matrix_view Phim = gsl_matrix_view_array(Phi, n, n);
    for (i = 0; i < n; i++)
      {
        double eval_i
           = gsl_vector_get (evalAp, i);
        gsl_vector_view evec_i
           = gsl_matrix_column (&Phim.matrix, i);

        //printf ("eigenvalue = %g\n", eval_i);
        //printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout,
                            &evec_i.vector, "%g");
      }
}


double find_2G(double *a)
{
    double pi15 = sqrt(M_PI*M_PI*M_PI);
    
    double H2[]={3.*pi15/(2.*sqrt(2.*a[0]))-3.*pi15/(sqrt(2.*a[0]))-M_PI/a[0], 
        6.*pi15*a[1]*a[1]/(pow(a[0]+a[1],2.5))-6.*pi15*a[1]/(pow(a[0]+a[1],1.5))-2.*M_PI/(a[0]+a[1]),
        6.*pi15*a[1]*a[1]/(pow(a[0]+a[1],2.5))-6.*pi15*a[1]/(pow(a[0]+a[1],1.5))-2.*M_PI/(a[0]+a[1]),
        3.*pi15/(2.*sqrt(2.*a[1]))-3.*pi15/(sqrt(2.*a[1]))-M_PI/a[1]
    };
    double S2[]={pi15/(pow(2.*a[0],1.5)),
                 pi15/(pow(a[0]+a[1],1.5)),
                 pi15/(pow(a[0]+a[1],1.5)),
                 pi15/(pow(2.*a[1],1.5))};
    
    double first_eig = gen_eigenvalues_nxn_gsl(H2, S2, 2); 
    printf("a1, a2 = %f, %f\teigenvalue = %f\n", a[0], a[1], first_eig);
    return first_eig;
}

double find_3G(double *a)   
{
    double pi15 = sqrt(M_PI*M_PI*M_PI);

    double H3[]={3.*pi15/(2.*sqrt(2.*a[0]))-3.*pi15/(sqrt(2.*a[0]))-M_PI/a[0], 
        6.*pi15*a[1]*a[1]/(pow(a[0]+a[1],2.5))-6.*pi15*a[1]/(pow(a[0]+a[1],1.5))-2.*M_PI/(a[0]+a[1]),
        6.*pi15*a[0]*a[0]/(pow(a[0]+a[2],2.5))-6.*pi15*a[0]/(pow(a[0]+a[2],1.5))-2.*M_PI/(a[0]+a[2]),
        6.*pi15*a[1]*a[1]/(pow(a[0]+a[1],2.5))-6.*pi15*a[1]/(pow(a[0]+a[1],1.5))-2.*M_PI/(a[0]+a[1]),
        3.*pi15/(2.*sqrt(2.*a[1]))-3.*pi15/(sqrt(2.*a[1]))-M_PI/a[1],
        6.*pi15*a[2]*a[2]/(pow(a[1]+a[2],2.5))-6.*pi15*a[2]/(pow(a[1]+a[2],1.5))-2.*M_PI/(a[1]+a[2]),
        6.*pi15*a[0]*a[0]/(pow(a[0]+a[2],2.5))-6.*pi15*a[0]/(pow(a[0]+a[2],1.5))-2.*M_PI/(a[0]+a[2]),
        6.*pi15*a[1]*a[1]/(pow(a[1]+a[2],2.5))-6.*pi15*a[1]/(pow(a[1]+a[2],1.5))-2.*M_PI/(a[1]+a[2]),
        3.*pi15/(2.*sqrt(2.*a[2]))-3.*pi15/(sqrt(2.*a[2]))-M_PI/a[2]
        };
        
    double S3[]={pi15/(pow(2.*a[0],1.5)),
                 pi15/(pow(a[0]+a[1],1.5)),
                 pi15/(pow(a[0]+a[2],1.5)),
                 pi15/(pow(a[0]+a[1],1.5)),
                 pi15/(pow(2.*a[1],1.5)),
                 pi15/(pow(a[1]+a[2],1.5)),
                 pi15/(pow(a[0]+a[2],1.5)),
                 pi15/(pow(a[1]+a[2],1.5)),
                 pi15/(pow(2.*a[2],1.5)),
        };
    
    double first_eig = gen_eigenvalues_nxn_gsl(H3, S3, 3);
    printf("a1, a2, a3 = %f, %f, %f\teigenvalue = %f\n", a[0], a[1], a[2], first_eig);
    return first_eig;
}



int main(int argc, char** argv)
{ 
    double pi15 = sqrt(M_PI*M_PI*M_PI);
    
    /*double a[] = {1.,5.,6.,
                  5.,1.,5.,
                  6.,5.,1.};
    double b[] = {2.,0.,1.,
                  0.,2.,0.,
                  1.,0.,2.};*/
    
    // Caso con singola gaussiana
    /*double Et=-0.5;
    double a1s=a1;
    
    while(abs(Et-E)>0.1){
        double H1[]={3.*pi15/(2.*sqrt(2.*a1s))-3.*pi15/(sqrt(2.*a1s))-M_PI/a1s};
        double S1[]={pi15/(pow(2.*a1s,1.5))};
            
        gen_eigenvalues_nxn_gsl(H1, S1,1); //generalized eigenvalue problem with gsl
    }*/
    // doppia gaussiana //NB PER RENDERE SIMMETRICA LA MATRICE USO IL FATTO CHE E' AUTOAGGIUNTA E AGISCE A SX QUANDO SERVE
    
    double a1d=a1; 
    double a2d=a2;
    double H2[]={3.*pi15/(2.*sqrt(2.*a1d))-3.*pi15/(sqrt(2.*a1d))-M_PI/a1d, 
        6.*pi15*a2d*a2d/(pow(a1d+a2d,2.5))-6.*pi15*a2d/(pow(a1d+a2d,1.5))-2.*M_PI/(a1d+a2d),
        6.*pi15*a2d*a2d/(pow(a1d+a2d,2.5))-6.*pi15*a2d/(pow(a1d+a2d,1.5))-2.*M_PI/(a1d+a2d),
        3.*pi15/(2.*sqrt(2.*a2d))-3.*pi15/(sqrt(2.*a2d))-M_PI/a2d
    };
    double S2[]={pi15/(pow(2.*a1d,1.5)),
                 pi15/(pow(a1d+a2d,1.5)),
                 pi15/(pow(a1d+a2d,1.5)),
                 pi15/(pow(2.*a2d,1.5))};
    
    gen_eigenvalues_nxn_gsl(H2, S2,2); 
    
    //Tripla gaussiana
    
    double a1t=a1; 
    double a2t=a2;
    double a3t=a3;
    double H3[]={3.*pi15/(2.*sqrt(2.*a1t))-3.*pi15/(sqrt(2.*a1t))-M_PI/a1t, 
        6.*pi15*a2t*a2t/(pow(a1t+a2t,2.5))-6.*pi15*a2t/(pow(a1t+a2t,1.5))-2.*M_PI/(a1t+a2t),
        6.*pi15*a1t*a1t/(pow(a1t+a3t,2.5))-6.*pi15*a1t/(pow(a1t+a3t,1.5))-2.*M_PI/(a1t+a3t),
        6.*pi15*a2t*a2t/(pow(a1t+a2t,2.5))-6.*pi15*a2t/(pow(a1t+a2t,1.5))-2.*M_PI/(a1t+a2t),
        3.*pi15/(2.*sqrt(2.*a2t))-3.*pi15/(sqrt(2.*a2t))-M_PI/a2t,
        6.*pi15*a3t*a3t/(pow(a2t+a3t,2.5))-6.*pi15*a3t/(pow(a2t+a3t,1.5))-2.*M_PI/(a2t+a3t),
        6.*pi15*a1t*a1t/(pow(a1t+a3t,2.5))-6.*pi15*a1t/(pow(a1t+a3t,1.5))-2.*M_PI/(a1t+a3t),
        6.*pi15*a2t*a2t/(pow(a2t+a3t,2.5))-6.*pi15*a2t/(pow(a2t+a3t,1.5))-2.*M_PI/(a2t+a3t),
        3.*pi15/(2.*sqrt(2.*a3t))-3.*pi15/(sqrt(2.*a3t))-M_PI/a3t
    };
    double S3[]={pi15/(pow(2.*a1t,1.5)),
                 pi15/(pow(a1t+a2t,1.5)),
                 pi15/(pow(a1t+a3t,1.5)),
                 pi15/(pow(a1t+a2t,1.5)),
                 pi15/(pow(2.*a2t,1.5)),
                 pi15/(pow(a2t+a3t,1.5)),
                 pi15/(pow(a1t+a3t,1.5)),
                 pi15/(pow(a2t+a3t,1.5)),
                 pi15/(pow(2.*a3t,1.5)),
    };
    
   
    gen_eigenvalues_nxn_gsl(H3, S3, 3);
    /*double aa[] = {1.,5.,6.,
                  5.,1.,5.,
                  6.,5.,1.};
    double bb[] = {2.,0.,1.,
                  0.,2.,0.,
                  1.,0.,2.};
    gen_eigenvalues_nxn_man(aa, bb, 3);
    //multiply3matrix(w, c, v, r, 3); function for three matrix product
    */
    
    double param[3];
    grad_descent(find_3G, 1e-6, 2.0, 3, param, 1);
    
    printf("\nmin = %f, %f, %f", param[0], param[1], param[2]);
    
    double param2d[3];
    grad_descent(find_3G, 1e-6, 1.0, 2, param2d, 1);
    printf("\nmin = %f, %f", param2d[0], param2d[1]);
    
    return 0;
}

/*void multiplytest( ){
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
    // now we do E=AD 
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &D.matrix,
                  0.0, &E.matrix);
    
    printf ("[ %g, %g, %g \n", e[0], e[1], e[2]);
    printf ("  %g, %g, %g \n", e[3], e[4], e[5]);
    printf ("  %g, %g, %g]\n", e[6], e[7], e[8]);
}*/
