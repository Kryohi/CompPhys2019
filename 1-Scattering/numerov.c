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
    printf ("Enter the number of gridpoints: ");
    scanf("%lf",&nmax);
    
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
    // creates array with the 1d eigenfunction
    double y[nmax];
    
    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++)
    {
            numerov(n, l, nmax, rmax, 0.1, y);
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

// Performs the whole algorithm and finds the spectrum up to the n-th level
// Returns 
double numerov(int nmax, int l, int xmax, double rmax, double Estep, double * y)
{
    double delta, prevdelta, xc_;
    int xc;
    int nfound = 0; // number of eigenvalues found;
    double h = rmax/nmax;
    double E = E0(h, rmax);
    double V_[xmax], k2[xmax];
    double yf[xmax], yb[xmax]; // ridurli alla lunghezza effettivamente usata?
    
    xc_ = secant(V, E, x1, x2, 1e-9, 1e-9) // finds the 0 of V(x) - E with the secant method. Should be faster, maybe.
    xc = (int) round(xc/h);
    
    for (int x==0; x<xmax; x++)
    {   V_[x] = V(x*h);
    }//if (V[x] >= E) xc = x;
    
    while (nfound <= nmax)
    {
        for (int x==0; x<xmax; x++)
            k2[x] = (E-V[x])/h2m - l*(l+1)/(x*h*x*h);
        

        delta = der3(yf,xc,h)/yf[xc] - der3(yb,xc,h)/yb[xc];
        
        // If there is a change in sign, start the finer search of the 0
        if (delta*prevdelta < 0)
        {
            E = zerosecant(double (*f)(double), double x1, double x2, double inf, double sup);
            
            nfound++;
        }
        prevdelta = delta;
        E += Estep;
    }
    
    return E;
}

/*
 * Differential equation iterative solution
 * aggiunge 2 punti dopo la barriera classica in xc, per calcolare derivata in quel punto
*/ 
void numerov_forward(int l, double E, double h, int xc, double * k2)
{
    y[0] = 0;
    y[1] = h^(l-1);

    for (int x=2; x<xc+2; x++)
        y[x] = (y[x-1]*(2 - 5*h*h*k2[x-1]/12) - y[x-2]*(1 + h*h*k2[x]/12)) / (1 + h*h*k2[x]/12); // controllare indici
    
    return y;
}


inline double V(double x)
{
    return x*x/2;
}


double E0(double h, double rmax)
{
    // we'll just assume that the minimum is at 0 ¯\_(ツ)_/¯
    return V(1e-9);
}



