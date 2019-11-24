/* hydrogen_gauss.c
   requires lapack 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lapacke.h>

#define ABS(a)     (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

main()
{
  /*
    Solution of the Hydrogen atom via expansion on a gaussian
    basis set and diagonalization (using linear algebra methods)
    NOTA BENE: fortran matrices are "simulated" by filling C vectors
               with the fortran ordering of matrix elements
  */

#define NALPX 100

  /* subroutine
  extern void diag (int , double *, double *, double *, double *); */
  /* variables */
  static double e2 = 2.0;
  static double pi = 3.14159265358979;
  int ia, ib, iab, n_alpha;
  double zeta, aa;
  double *alpha, *e, *s, *h, *v; 
  double r, f, prob, q, dr;
  int i, j, nrx;
  char filin[80];
  FILE *in, *out;

  /* Input data */

  fprintf(stdout, " Atomic charge >> ");
  scanf("%lf", &zeta);
  fprintf(stdout, " Parameters of the Gaussians from file >> ");
  scanf("%80s", filin);
  in = fopen (filin,"r");
  if ( in == NULL ) {
    fprintf(stderr, "file %s not found\n",filin);
    exit (1);
  } else {
    fprintf(stdout, "\n");
  }
  fscanf(in,"%d", &n_alpha);
  printf("nalpha = %d\n",n_alpha);
  if ( n_alpha > NALPX)  {
    fprintf(stderr, "n_alpha (%d) > nalpx (%d)\n",n_alpha, NALPX);
    exit (1);
  } else {
    alpha = (double *) malloc ( n_alpha * sizeof (double) );
  }
  /* Read parameters of the gaussians */

  for ( ia = 0 ; ia < n_alpha ; ++ia ) {
    fscanf(in,"%lf", &alpha[ia]);
  }
    printf("alpha0 = %f\t alpha1 = %f\n",alpha[0],alpha[1]);
  fclose(in);

  e = (double *) malloc ( n_alpha * sizeof (double) );
  h = (double *) malloc ( n_alpha * n_alpha * sizeof (double) );
  v = (double *) malloc ( n_alpha * n_alpha * sizeof (double) );
  s = (double *) malloc ( n_alpha * n_alpha * sizeof (double) );

  /* ***************************************************************
   S States
  ****************************************************************** */

  /* fprintf (stdout, "Hydrogen-like problem: 1 electron, S states\n"); */

  /* Assign values of overlap integrals and of matrix elements 
     of the hamiltonian on the S basis: */
  iab = 0;
  for ( ib = 0 ; ib < n_alpha ; ++ib ) {
    for ( ia = 0 ; ia < n_alpha ; ++ia ) {
      aa = alpha[ia] + alpha[ib];
      s[iab] = pow( (pi/aa), 1.5 );
      h[iab] = s[iab]*6.0*alpha[ia]*alpha[ib]/aa - e2*zeta*2.*pi/aa;
      iab++;
    }
  }

  /* Solution [expansion coefficients are stored into v(j,i)
     j=basis function index, i= eigenvalue index] */

  diag ( n_alpha, h, s, e, v ) ;

  fprintf (stdout, "Lowest eigenvalues, S states: %f %f %f\n", e[0],e[1],e[2] );

  /* Write to output file the coefficients of the lowest-energy state: */

  out = fopen ("s-coeff.out","w");

  for ( j = 0 ; j < n_alpha ; ++j ) {
    fprintf (out,"%d %f %f\n", j, alpha[j], v[j] ) ;
  }
  fclose(out);

  /*   Write to output file the lowest-energy state: */
  
  out = fopen ("s-wfc.out","w");
  dr  = 0.1/zeta;
  nrx = 200;
  q   = 0.0;
  for ( i = 0 ; i <= nrx; ++i ) {
    r = dr*i;
    f = 0.0;
    for ( j = 0; j < n_alpha; ++j ) {
      f = f + exp(-alpha[j]*r*r)*v[j];
    }
    prob = 4.0 * pi * pow ( (r*f), 2);
    q   = q + prob*dr ;
    fprintf (out,"%f %f %f %f\n", r, f, r*f, prob) ;
  }
  /* verify normalization (if desired):
     write (*,*) q */

  /* ***************************************************************
   P States
  ****************************************************************** */

 /* fprintf (stdout, "Hydrogen-like problem: 1 electron, P states\n"); */

  /* Assign values of overlap integrals and of matrix elements 
     of the hamiltonian on the P basis: */
  iab = 0;
  for ( ib = 0 ; ib < n_alpha ; ++ib ) {
    for ( ia = 0 ; ia < n_alpha ; ++ia ) {
      aa = alpha[ia] + alpha[ib];
      s[iab] = 0.5 / aa * pow( (pi/aa), 1.5 );
      h[iab] = s[iab]*10.0*alpha[ia]*alpha[ib]/aa - 
                  e2*zeta*2.*pi/aa/aa/3.0;
      iab++;
    }
  }

  /* Solution [expansion coefficients are stored into v(j,i)
     j=basis function index, i= eigenvalue index] */

  diag ( n_alpha, h, s, e, v ) ;

  fprintf (stdout, "Lowest eigenvalues, P states: %f %f %f\n", e[0],e[1],e[2] );

  /* Write to output file the coefficients of the lowest-energy state: */

  out = fopen ("p-coeff.out","w");

  for ( j = 0 ; j < n_alpha ; ++j ) {
    fprintf (out,"%d %f %f\n", j, zeta/sqrt(alpha[j]), v[j] ) ;
  }
  fclose(out);

  /*   Write to output file the lowest-energy state: */
  
  out = fopen ("p-wfc.out","w");
  dr  = 0.1/zeta;
  nrx = 200;
  q   = 0.0;
  for ( i = 0 ; i <= nrx; ++i ) {
    r = dr*i;
    f = 0.0;
    for ( j = 0; j < n_alpha; ++j ) {
      f = f + r * exp(-alpha[j]*r*r)*v[j];
    }
    prob = 4.0/3.0 * pi * pow ( (r*f), 2);
    q   = q + prob*dr ;
    fprintf (out,"%f %f %f %f\n", r, f, r*f, prob) ;
  }
  /* verify normalization (if desired):
     write (*,*) q */
  fclose(out);
  free(s); free(v); free(h); free(e); free (alpha);
}





/* subroutine diag */
diag (int n, double *h, double *s, double *e, double *v)
{
  /*    Finds eigenvalues and eigenvectors of the generalized problem
	Hv=eSv, where H=hermitian matrix, S=overlap matrix */

  /* On input: n = dimension of the matrix to be diagonalized
               h = matrix to be diagonalized
               s = overlap matrix
     On output:
               e = eigenvalues
               v = eigenvectors
               s and h are unchanged */

  /* LOCAL variables */
  int lwork, i, j, k, nn, ij, info;
  /* lwork = dimension of workspace for lapack routine  */
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, *aux, *uno;

  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  /* Copy S into an auxiliary matrix (dsyev destroys the matrix) */
  aux = (double *) malloc( n * n * sizeof (double));
  for ( ij = 0; ij < n*n; ++ij ) {
    aux[ij] = s[ij];
  }

  /*  Diagonalize S  */

  dsyev_ ( V, U, &n, aux, &n, e, work, &lwork, &info ) ;

  if ( info !=0 ) {
    fprintf (stderr, "S-matrix diagonalization failed\n");
    exit (1);
  }

  /*    Keep only linearly independent combinations
	(within a given threshold)  */

  nn = 0;
  for ( i = 0; i < n; ++i ) {
    /*  i runs on all eigenstates
       nn runs on eigenstate with nonzero eigenvalue */
    if ( e[i] > small) {
      for ( j = 0; j < n; ++j ) {
	aux[j+nn*n] = aux[j+i*n] / sqrt(e[i]);
      }
      ++nn;
    }
  }

  if ( nn < n ) { 
     fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
  }
  /*       Trasform H using the "aux" matrix
           V(i,j) = \sum_{k=1}^{n} H(i,k) aux(k,j),  i=1,n, j=1,nn
   */
  
  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
      }
    }
  }
   /* h1(i,j) = \sum_{k=1}^{n} aux(k,i) v(k,j),  i=1,nn, j=1,nn
      H' = transpose(aux)*H*aux
   */
  uno = (double *) malloc( nn * nn * sizeof (double));
  for ( i = 0; i < nn; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      uno[i+j*nn] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
      }
    }
  }
  
  /*  Diagonalize transformed H  */
  
  dsyev_ ("V", "U", &nn, uno, &nn, e, work, &lwork, &info );

  if ( info !=0 ) {
    fprintf (stderr, "H-matrix diagonalization failed\n");
    exit (1);
  }

  /*  Back-transform eigenvectors  */

  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < nn; ++k ) {
	v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
      }
    }
  }
  free(uno); free(aux); free(work);
  return;
}
