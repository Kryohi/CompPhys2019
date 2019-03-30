#include "numerov.c"
#include "misccose.c"

double V_lj(double x, double epsilon, double sigma);
double V_fastlj(double dr);

#define lmax 6 // 6
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
    double Emax = 0.59322; // 3.5/5.9
    double m = 1.66054e-27/(1/83.798 + 1/1.00794); // in kg
    double h2m = hbar*hbar/(2*m); // in un ibrido bastardizzato SI/eV
    h2m = h2m/(epsilon*sigma*sigma); // ridotta
    double rmax = 7; // 7*sigma o 7 e basta in unità naturali ?
    int xmax = 5000; //(int)round(rmax/h);
    double h = rmax/xmax;
    
    double k, a, sigma_tot;
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
    double b = pow(4/(25*h2m), 0.1);
    printf("b = %f\t", b); // ~1.1635
    for (int x=1; x<bc.length; x++)
        bc.data[x] = exp(-pow(b/(x*h), 5));  // mi son dimenticato cos'era b, il valore di V quando passa da 0?
    printf("bc[%d] = %f\n", xmin-1, bc.data[xmin-1]);
    
    /*    ROBA INUTILE che vorrei calcolare solo per curiosità 
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
    
    
    /* Scattering, TODO */
    
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
        
    // Bessel functions in r1
    J1[0] = sin(r1)/r1;
    J1[1] = (J1[0] - cos(r1))/r1;
    N1[0] = -cos(r1)/r1;
    N1[1] = (N1[0] - sin(r1))/r1;
    fast_bessel(k*r1, lmax, J1);
    fast_bessel(k*r1, lmax, N1);
    
    // Bessel functions in r2
    J2[0] = sin(r2)/r2;
    J2[1] = (J2[0] - cos(r2))/r2;
    N2[0] = -cos(r2)/r2;
    N2[1] = (N2[0] - sin(r2))/r2;
    fast_bessel(k*r2, lmax, J2);
    fast_bessel(k*r2, lmax, N2);
    
    k = sqrt(E/h2m);
    a = y[x1]*r2 / (y[x2]*r1);

    for (int l=0; l<=lmax; l++)
        delta[l] = arctan( (c*) / () )
    
    sigma_tot = (4*MATH_PI / (k*k)*sum((2*l+1)*sin(delta_l)*sin(delta_l));


    free(bc.data);
    
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

