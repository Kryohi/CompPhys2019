#include "numerov.c"
#include "misccose.c"

#define lmax 6
#define rlow 
#define h2m 0.5

/* TODO 
 * modificare bc[] in modo che possa avere lungh arbitraria (in teoria fatto)
 * oppure, far partire numerov da nuovo parametro xlow
 * precalcolare funzione del punto 5 fino a rlow e passarla a Numerov / aggiungerla a autofunzione finale
 * 
 */

int main(int argc, char** argv)
{
    
    double rmax = 7*sigma;
    int xmax = (int)round(rmax/h);
    
    double c, sigma_tot;
    double delta[lmax];
    double J1[lmax], double N1[lmax], double J2[lmax], double N2[lmax];
    
    // creates data folder and files to save data
    make_directory("Data_scattering");
    chdir("Data_scattering");
    
    // Boundary conditions near zero, see point 5
    int xmin = (int)round(rlow/h);   // rlow in the notes, point where V equals 100
    dArray bc;
    bc.length = xmin;
    bc.data = calloc(xmin, sizeof(double));
    bc.data[0] = 0;
    for (int x=1; x<bc.length; x++)
        bc.data[x] = exp(-pow(x*h/b, -5));
    
    // creates array of structs with the results, one set for each value of l
    struct Spectrum spectra[lmax+1];
    for (int l=0; l<=lmax; l++) {
        spectra[l].xmax = xmax;
        spectra[l].nmax = nmax;
        spectra[l].EE = calloc(nmax, sizeof(double));
        spectra[l].eigfuns = calloc(xmax*nmax, sizeof(double));
    }
    
    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++) {
        //bc.data[1] = pow(rmax/xmax,l);
        spectra[l] = numerov(nmax, l, xmax, rmax, 0.13, true, bc, V_ho);
    }
    
    
    double r1 = 5.5*sigma;
    double r2 = 6.0*sigma;
    int x1 = (int)round(r1/h);
    int x2 = (int)round(r2/h);
    
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



inline double V_lj(double x, double epsilon, double sigma)
{
    return 4*epsilon*(pow(sigma/x,12)-pow(sigma/x,6));
}

inline double V_fastlj(double dr2)
{
     double dr6 = dr2*dr2*dr2;
     return 4*(1.0/(dr6*dr6) - 1.0/dr6);
}

