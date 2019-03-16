// TODO
/*
 * subito con matching all'indietro?
 * come strutturare output
 * cosa serve cambiare per scattering
 * 
 */

#include "numerov.h"
#include "matematicose.c"
#include "misccose.c"

// hbar/2m, to change to LJ natural units at point 4
#define h2m 0.5

int main(int argc, char** argv)
{
    int nmax, lmax, xmax;
    double rmax;
    
    // asks user for grid parameters and quantum numbers
    printf("Enter the maximum quantum number n:");
    scanf("%d",&nmax);
    printf("Enter the maximum quantum number l:");
    scanf("%d",&lmax);
    printf("Enter rmax: ");
    scanf("%lf",&rmax);
    printf("Enter the number of gridpoints: ");
    scanf("%d",&xmax);
    
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
    // creates struct with the results
    //double y[nmax];
    
    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++)
    {
        numerov(nmax, l, xmax, rmax, 0.1);
    }
    
    
    // save to a csv file    
    char filename[64];
    snprintf(filename, 64, "./eigen_%d.csv", nmax);
    FILE * eigen;
    eigen = fopen(filename, "w");
    if (eigen == NULL)
        perror("error while writing on eigen.csv");
    
    return 0;
}

// Performs the whole algorithm and finds the spectrum up to the n-th level
// Returns 
double numerov(int nmax, int l, int xmax, double rmax, double Estep)
{
    double delta, delta1, delta2, prevdelta=0;
    double xc_float;
    int xc; // point of classical barrier
    int nfound = 0; // number of eigenvalues found;
    int niter = 0;
    double h = rmax/nmax;
    double E=E0(h, rmax), E1, E2, E_;
    double V_[xmax], k2[xmax];
    double yf[xmax], yb[xmax]; // ridurli alla lunghezza effettivamente usata?
    
    for (int x=0; x<xmax; x++)
        V_[x] = V(x*h);
        
    // Boundary conditions (TODO)
    yf[0] = 0;
    yf[1] = pow(h, l-1);
    yb[xmax-1] = 0;
    yb[xmax-2] = 1e-10;
    
    
    while (nfound <= nmax)
    {
        xc_float = secant(V, E, 0., rmax, h/4, h/4); // finds the 0 of V(x)-E with the secant method.
        xc = (int) round(xc_float/h);
        for (int x=0; x<xmax; x++)  k2[x] = (E-V_[x])/h2m - l*(l+1)/(x*h*x*h);
        numerov_forward(h, xc, k2, yf);
        numerov_backward(h, xc, xmax, k2, yb);
        delta = der3(yf,xc,h)/yf[xc] - der3(yb,xc,h)/yb[xc];
        
        // If there is a change in sign, start the finer search of the 0 with the secant method
        if (delta*prevdelta < 0)
        {
            E1 = E - Estep; // store the previous value of E, before the change of sign
            E2 = E;
            delta1 = prevdelta;
            delta2 = delta;
            
            while (fabs(delta2) < 1e-4)
            {
                E_ = E2;  // E precedente
                E2 = E2 - delta2*(E2-E1) / (delta2-delta1);
                E1 = E_;
                delta1 = delta2;
                
                xc_float = secant(V, E2, 0., rmax, h/4, h/4); // finds the 0 of V(x)-E with the secant method.
                xc = (int) round(xc_float/h);
                for (int x=0; x<xmax; x++)  k2[x] = (E2-V_[x])/h2m - l*(l+1)/(x*h*x*h);
                numerov_forward(h, xc, k2, yf);
                numerov_backward(h, xc, xmax, k2, yb);
                delta2 = der3(yf,xc,h)/yf[xc] - der3(yb,xc,h)/yb[xc];
                
                niter++;
            }
            printf("E%d = %f\n", nfound, E2);
            printf("Found in %d iterations\n\n", niter);
            niter = 0;
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
void numerov_forward(double h, int xc, const double * k2, double * y)
{
    for (int x=2; x<xc+2; x++)
        y[x] = (y[x-1]*(2 - 5*h*h*k2[x-1]/12) - y[x-2]*(1 + h*h*k2[x]/12)) / (1 + h*h*k2[x]/12); // controllare indici
}

void numerov_backward(double h, int xc, int xmax, const double * k2, double * y)
{
    for (int x=xmax-3; x<xc-2; x--)
        y[x] = (y[x+1]*(2 - 5*h*h*k2[x+1]/12) - y[x+2]*(1 + h*h*k2[x]/12)) / (1 + h*h*k2[x]/12); // controllare indici
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



