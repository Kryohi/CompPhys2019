#include "numerov.c"
#include "misccose.c"

#define lmax 6

/* TODO 
 * modificare potenziale in modo da spostare barriera + vicina a 0? difficile perché in numerov entra la x
 * oppure, meglio, modificare bc[] in modo che possa avere lungh arbitraria
 * oppure, far partire numerov da nuovo parametro xlow
 * precalcolare funzione del punto 5 fino a rlow e passarla a Numerov / aggiungerla a autofunzione finale
 * 
 */

int main(int argc, char** argv)
{
    
    double rmax = 7*sigma;
    
    double c, sigma_tot;
    double delta[lmax];
    double J1[lmax], double N1[lmax], double J2[lmax], double N2[lmax];
    
    
    
    c = y[x1]*x2 / (y[x2]*x1*h);
    fast_bessel(k*r1, J1);
    for (int l=0; l<=lmax; l++)
        delta[l] = arctan( (c*) / () )
    
    sigma_tot = (4*MATH_PI / (k*k)*sum((2*l+1)*sin(delta_l)*sin(delta_l));
    
    return 0;
}

double phase_shift()
{
    
    
}
