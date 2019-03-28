#include "numerov.c"
#include "misccose.c"

#define lmax 6
#define rlow 0.4 //TODO
#define hbar 1.0545718e-34
#define hbar_eV 4.13566766225e-15

/* TODO 
 * 
 * 
 */

int main(int argc, char** argv)
{
    double sigma = 3.18e-10;
    double epsilon = 5.9e-3 * 1.602176565e-19;
    double Emax = 0.59322; // 3.5/5.9
    double m = 1.66054e-27/(1/83.798 + 1/1.00794); // prob da mettere in grammi
    double h2m = hbar*hbar/(2*m); // in un ibrido bastardizzato SI/eV
    h2m = h2m/(epsilon*sigma*sigma); // ridotta
    double rmax = 7; // 7*sigma o 7 e basta in unità naturali ?
    int xmax = 10000; //(int)round(rmax/h);
    double h = rmax/xmax;
    
    double c, sigma_tot;
    double delta[lmax];
    double J1[lmax]; double N1[lmax]; double J2[lmax]; double N2[lmax];
    
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
        bc.data[x] = exp(-pow(x*h/b, -5));  // mi son dimenticato cos'era b, il valore di V quando passa da 0?

    
    /*    ROBA INUTILE che vorrei calcolare solo per curiosità     */
    
    // creates array of structs with the results, one set for each value of l
    struct Spectrum spectra[lmax+1];
    for (int l=0; l<=lmax; l++) {
        spectra[l].xmax = xmax;
        spectra[l].nmax = 7;
        spectra[l].EE = calloc(nmax, sizeof(double));
        spectra[l].eigfuns = calloc(xmax*nmax, sizeof(double));
    }
    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++) {
        spectra[l] = numerov(7, l, xmax, rmax, h2m, 0.13, true, bc, V_fastlj);
    }
    // saves to two csv files all the energy levels and the eigenfunctions
    save2csv(spectra, lmax, nmax, xmax);
    
    
    
    /* Scattering, TODO
     * 
    
    // precalculation of the known terms
    for (int x=0; x<xmax; x++)
    {
        V_[x] = V_fastlj(x*h);
        centrifugal[x] = l*(l+1)/(x*h*x*h+1e-14);
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
    */

    free(bc.data);
    for (int l=0; l<=lmax; l++)
        free(spectra[l].eigfuns);
    return 0;
}

double phase_shift()
{
    
    
}



inline double V_lj(double x, double epsilon, double sigma)
{
    return 4*epsilon*(pow(sigma/x,12)-pow(sigma/x,6));
}

inline double V_fastlj(double dr)
{
     double dr6 = dr*dr*dr*dr*dr*dr;
     return 4*(1.0/(dr6*dr6) - 1.0/dr6);
}

