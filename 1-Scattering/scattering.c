#include "numerov.c"
#include "misccose.c"

double V_lj(double x, double epsilon, double sigma);
double V_fastlj(double dr);

#define lmax 9 // 6
#define rlow 0.63 // where V = 10 0.923, V=100 0.752, V=1000 0.63
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
    double E, Emax = 0.59322; // 3.5/5.9
    double m = 1.66054e-27/(1/83.798 + 1/1.00794); // in kg
    double h2m = hbar*hbar/(2*m); // in un ibrido bastardizzato SI/eV
    h2m = h2m/(epsilon*sigma*sigma); // ridotta
    double rmax = 10.01; // Melius est abundare quam deficere
    int xmax = 16000; // vedi sopra
    double h = rmax/xmax;
    int xmin = (int)round(rlow/h);   // point where V equals 1000 (chosen by looking at when the bound states become uncorrelated to this)
    printf("\nxmin is %d\n", xmin);
    double y[xmax], k2[xmax], centrifugal[xmax], V_[xmax];
    
    double k, a, sigma_tot = 0;
    double delta[lmax+1];
    double J1[lmax+1]; double N1[lmax+1]; double J2[lmax+1]; double N2[lmax+1];
    
    // creates data folder and files to save data
    make_directory("Data_scattering");
    chdir("Data_scattering");
    
    // Files for output
    FILE * wf;
    FILE * ps;
    FILE * tcs;
    char filename[64];
    
    
    
    /*    ROBA INUTILE che calcolo solo per curiosità 
     dArray bc;
    bc.length = xmin;
    bc.data = calloc(xmin, sizeof(double));
    bc.data[0] = 0;
    double b = pow(4/(25*h2m), 0.1);
    printf("b = %f\t", b); // ~1.1635
    
    for (int x=1; x<bc.length; x++)
        bc.data[x] = exp(-pow(b/(x*h), 5));  // mi son dimenticato cos'era b, il valore di V quando passa da 0?
        
    printf("bc[%d] = %f\n", xmin-1, bc.data[xmin-1]);
    int nmax = 1; // esiste solo uno stato legato (?)
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
        spectra[l] = numerov(nmax, l, xmax, rmax, h2m, 0.01, true, bc, V_fastlj);
    }
    // saves to two csv files all the energy levels and the eigenfunctions
    save2csv(spectra, lmax, nmax, xmax);
    */
    
    
    /*   Scattering   */
    
    double r1 = 5.0;//*sigma;
    double r2 = r1+.5;//*sigma;
    double kr1, kr2;
    int x1 = (int)round(r1/h);
    int x2 = (int)round(r2/h);
    
    // Energy fixed for point 6
    E = 0.3;

    
    // precalculation of the known, constant terms
    
        k = sqrt(E/h2m);
        kr1 = k*r1;
        kr2 = k*r2;
        
        // Bessel functions in r1
        J1[0] = sin(kr1)/kr1;
        J1[1] = (J1[0] - cos(kr1))/kr1;
        N1[0] = -cos(kr1)/kr1;
        N1[1] = (N1[0] - sin(kr1))/kr1;
        fast_bessel(kr1, lmax, J1);
        fast_bessel(kr1, lmax, N1);
    
        // Bessel functions in r2
        J2[0] = sin(kr2)/kr2;
        J2[1] = (J2[0] - cos(kr2))/kr2;
        N2[0] = -cos(r2)/r2;
        N2[1] = (N2[0] - sin(kr2))/kr2;
        fast_bessel(kr2, lmax, J2);
        fast_bessel(kr2, lmax, N2);
        
        for (int l=0; l<=lmax; l++)  {
            printf("J1[%d] = %f\n", l, J1[l]);
            printf("N1[%d] = %f\n", l, N1[l]);
            printf("J2[%d] = %f\n", l, J1[l]);
            printf("N2[%d] = %f\n", l, N2[l]);
        }printf("\n");

    for (int x=0; x<xmax; x++)  {
        V_[x] = V_fastlj(x*h+1e-15);
        centrifugal[x] = 1/(x*h*x*h+1e-15);
    }
    
    // Boundary conditions near zero, see point 5
    y[0] = 0;
    double b = pow(4/(25*h2m), 0.1); // ~1.1635, da ricontrollare?
    printf("b = %f\n", b);
    for (int x=1; x<xmin; x++)
        y[x] = exp(-pow(b/(x*h), 5));
    
    
    for (int l=0; l<=lmax; l++)
    {
        // Numerov solution
        for (int x=0; x<xmax; x++)
            k2[x] = (E-V_[x])/h2m - l*(l+1)*centrifugal[x];
    
        numerov_forward(h*h, xmax, xmin, k2, y);
    
        a = y[x1]*r2 / (y[x2]*r1);
        delta[l] = atan2( (a*J2[l] - J1[l]), (a*N2[l] - N1[l]) );   // phase shift
        sigma_tot += 4*M_PI * (2*l+1) * sin(delta[l])*sin(delta[l]) / (k*k);
        
        snprintf(filename, 64, "./wavefunction_l%d.csv", l);
        wf = fopen(filename, "w");
        fprintf(wf, "Veff, y\n");
        for (int x=0; x<xmax; x++)
            fprintf(wf, "%f, %0.12f\n", V_[x] + h2m*l*(l+1)*centrifugal[x], y[x]);
    }
    
    for (int l=0; l<=lmax; l++)
        printf("phase shift l=%d: %f\t sin^2 = %f\n", l, delta[l], sin(delta[l])*sin(delta[l]));
    
    printf("scattering amplitude : %f\n", sigma_tot);
    
    snprintf(filename, 64, "./phaseshifts_rlow%0.4f_rmax%0.4f.csv", rlow, r2);
    ps = fopen(filename, "w");
    fprintf(ps, "ph_sh, contribution, sigma_tot\n");
    for (int l=0; l<=lmax; l++)
        fprintf(ps, "%0.9f, %0.9f, %f\n", delta[l], (2*l+1) * sin(delta[l])*sin(delta[l]), sigma_tot);
    
    
    return 0;
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

