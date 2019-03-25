#include "numerov.c"

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

    // Iterates the Numerov algorithm for different quantum numbers
    for (int l=0; l<=lmax; l++)
        spectra[l] = numerov(nmax, l, xmax, rmax, 0.13, true, V_ho);
    
    
    // saves to two csv files all the energy levels and the eigenfunctions
    save2csv(spectra, lmax, nmax, xmax);
    
    return 0;
}
