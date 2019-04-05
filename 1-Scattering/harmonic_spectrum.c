#include "numerov.c"

#define h2m 0.5


int main(int argc, char** argv)
{
    int nmax, lmax, xmax;
    double rmax;
    
    if (argc == 5)
    {
        lmax = (int)strtol(argv[1], NULL, 10);
        nmax = (int)strtol(argv[2], NULL, 10);
        rmax = (float)strtol(argv[3], NULL, 10);
        xmax = (int)strtol(argv[4], NULL, 10);
    }
    else    {
        // asks user for grid parameters and quantum numbers
        printf("Enter the maximum quantum number l: ");
        scanf("%d",&lmax);
        printf("Enter the maximum quantum number n: ");
        scanf("%d",&nmax);
        printf("Enter rmax: ");
        scanf("%lf",&rmax);
        printf("Enter the number of gridpoints: ");
        scanf("%d",&xmax);
    }
    
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
    // creates array of structs with the results, one set for each value of l
    struct Spectrum spectra[lmax+1];
    for (int l=0; l<=lmax; l++) {
        spectra[l].xmax = xmax;
        spectra[l].nmax = nmax;
        spectra[l].EE = calloc(nmax, sizeof(double));
        spectra[l].eigfuns = calloc(xmax*nmax, sizeof(double));
    }
    
    // Boundary conditions near zero, they depend on V
    dArray bc;
    bc.length = 2;
    bc.data = calloc(2, sizeof(double));
    bc.data[0] = 0;
    
    
    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++) {
        bc.data[1] = pow(rmax/xmax,l);
        spectra[l] = numerov(nmax, l, xmax, rmax, h2m, 0.13, true, bc, V_ho);
    }
    
    
    // saves to two csv files all the energy levels and the eigenfunctions
    save2csv(spectra, lmax, nmax, xmax);
    
    free(bc.data);
    for (int l=0; l<=lmax; l++)
        free(spectra[l].eigfuns);
        
    return 0;
}



inline double V_ho(double x)
{
    int ho1D = 0;
    if (ho1D==1) {
        return (x-10.)*(x-10.)/2;
    }
    else {
        return (x)*(x)/2;
    }
}

