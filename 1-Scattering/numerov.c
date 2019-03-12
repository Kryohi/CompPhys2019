// TODO
/*
 * algoritmo numerov
 * subito con matching all'indietro?
 * come strutturare output
 * cosa serve cambiare per scattering
 * 
 */

#include "numerov.h"
#include "matematicose.c"

// hbar/2m, to change to LJ natural units at point 4
#define h2m 0.5

int main(int argc, char** argv)
{
    int nmax, lmax;
    double rmax, h, r, unorm, uu;
    
    // asks user for grid parameters
    printf ("Enter rmax: ");
    scanf("%lf",&rmax);
    printf ("Enter h: ");
    scanf("%lf",&h);
    
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
    // creates array with the 1d eigenfunction
    double x [ceil(rmax/h)];
    
    // Iterates the Numerov algorithm for different quantum numbers
    for (int n=0; n<=nmax; n++)
    {
        for (int l=0; l<=lmax; l++)
        {
            numerov(n, l, h, rmax, x);
        }
    }
    
    
    
    // save to a csv file    
    char filename[64];
    snprintf(filename, 64, "./eigen_N%d_M%d_r%0.4f_T%0.2f.csv", N, M, rho, T);
    FILE * eigen;
    eigen = fopen(filename, "w");
    if (eigen == NULL)
        perror("error while writing on eigen.csv");
    
    return 0;
}

void numerov(int n, int l, double h, double rmax, double * y)
{
    k2 = (1/h2m)*(E-V(x)) - l*(l+1)/x^2;
    
    for (int x=2; x <= ceil(rmax/h); x++)
    {
        y[x] = (y[x-1]*(2 - 5*h*h*k2[x-1]/12) - y[x-2]*(1 + h*h*k2[x]/12)) / (1 + h*h*k2[x]/12); // controllare indici
    }
    
}


inline double V(const double x)
{
    return x*x/2;
}

