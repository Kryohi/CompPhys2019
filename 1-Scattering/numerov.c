// TODO
/*
 * fix inversione di segno appena dopo autovalori
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
        spectra[l] = numerov(nmax, l, xmax, rmax, 0.13, V_ho);
    
    
    // saves to two csv files all the energy levels and the eigenfunctions
    save2csv(spectra, lmax, nmax, xmax);
    
    return 0;
}

// Performs the whole algorithm and finds the spectrum up to the n-th level
// Returns a Spectrum struct
Spectrum numerov(int nmax, int l, int xmax, double rmax, double Estep, double (*V)(double))
{
    Spectrum sp;
    sp.EE = calloc(nmax, sizeof(double));
    sp.eigfuns = calloc(xmax*nmax, sizeof(double));
    double h = rmax/xmax;
    double xc_float;
    int xc; // point corresponding to the classical barrier
    int nfound = 0; // number of eigenvalues found;
    double E=E0(h, rmax, V), E1, E2, E_;
    double delta, delta1, delta2, prevdelta=0;
    double V_[xmax], centrifugal[xmax], k2[xmax];
    double yf[xmax], yb[xmax];
    printf("\nFinding the spectrum for l=%d...\n", l);
    
    for (int x=0; x<xmax; x++)
    {
        V_[x] = V(x*h);
        centrifugal[x] = l*(l+1)/(x*h*x*h+1e-14);
    }
    
    // Boundary conditions near 0
    yf[0] = 0;
    yf[1] = pow(h,l);
    
    while (nfound < nmax)
    {
        //xc_float = secant(V, E, 1e-6, rmax, -h/4, h/4); // finds the 0 of V(x)-E with the secant method.
        xc_float = sqrt(2*E); // !!!! vale solo per armonico !!!!
        xc = (int) round(xc_float/h);

        for (int x=0; x<xmax; x++)  
            k2[x] = (E-V_[x])/h2m - centrifugal[x];
        
        // Boundary conditions at rmax //TODO controllare se ci sono da imporre segni alternanti
        yb[xmax-1] = exp(-sqrt(fabs(E)/h2m)*xmax*h);
        yb[xmax-2] = exp(-sqrt(fabs(E)/h2m)*(xmax-1)*h);
        
        numerov_forward(h*h, xc, k2, yf);
        numerov_backward(h*h, xc, xmax, k2, yb);
        //printf("yf[xc-1] = %f,  yf[xc] = %f,  yf[xc+1]=%f\n", yf[xc-1], yf[xc], yf[xc+1]);
        //printf("yb[xc-1] = %f,  yb[xc] = %f,  yb[xc+1]=%f\n", yb[xc-1], yb[xc], yb[xc+1]);
        delta = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc];
        //printf("derforward = %f,  derback = %f\n", der3(yf,xc,h)/yf[xc], der3(yb,xc,h)/yb[xc]);
        printf("E = %f\tDelta rough = %f\n", E, delta);
        
        // If there is a change in sign, start the finer search of the 0 of delta(E) with the secant method
        if (delta*prevdelta < 0)
        {
            printf("\nFound a point of inversion!\n");
            E1 = E - Estep; // store the previous value of E, before the change of sign
            E2 = E;
            delta1 = prevdelta;
            delta2 = delta;
            
            while (fabs(delta2) > 1e-6)
            {
                // Secant method
                E_ = E2;  // E precedente
                E2 = E2 - delta2*(E2-E1) / (delta2-delta1);
                E1 = E_;
                delta1 = delta2;
                printf("E = %f\tdelta2 = %f\n", E2, delta2);
                
                // Calculation of the next delta2 value
                xc_float = secant(V, E2, 0., rmax, -h/4, h/4); // finds the 0 of V(x)-E with the secant method.
                xc_float = sqrt(2*E); // !!!! vale solo per armonico !!!!
                xc = (int) round(xc_float/h);
                for (int x=0; x<xmax; x++)
                    k2[x] = (E2-V_[x])/h2m - centrifugal[x];
                
                yb[xmax-1] = exp(-sqrt(fabs(E)/h2m)*xmax*h);
                yb[xmax-2] = exp(-sqrt(fabs(E)/h2m)*(xmax-1)*h);
                numerov_forward(h*h, xc, k2, yf);
                numerov_backward(h*h, xc, xmax, k2, yb);
                delta2 = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc];
                printf("Delta fine = %f\n", delta2);
            }
            
            printf("E%d = %.9f\n\n", nfound, E2);
            sp.EE[nfound] = E2;
            for (int x=0; x<xc; x++)
                sp.eigfuns[xmax*nfound + x] = yf[x];
            for (int x=xc; x<xmax; x++)
                sp.eigfuns[xmax*nfound + x] = yb[x] * yf[xc]/yb[xc]; //impose continuity
            nfound++;
        }
        
        prevdelta = delta;
        E += Estep;
        
        if (E>V_[xmax-1])
            perror("Could not find enough solutions with the parameters provided\n");
    }
    
    return sp;
}

/*
 * Differential equation iterative solution
 * aggiunge 2 punti dopo la barriera classica in xc, per calcolare derivata in quel punto
*/ 
void numerov_forward(double hh, int xc, const double * k2, double * y) 
{
    double c0 = hh*k2[2]/12;
    double c_1 = hh*k2[1]/12;
    double c_2 = hh*k2[0]/12;
    for (int x=2; x<xc+3; x++)  {
        y[x] = (y[x-1]*(2-10*c_1) - y[x-2]*(1+c_2)) / (1+c0);
        c_2 = c_1;
        c_1 = c0;
        c0 = hh*k2[x+1]/12;
    }
}

void numerov_backward(double hh, int xc, int xmax, const double * k2, double * y)
{
    double c0 = hh*k2[xmax-3]/12;
    double c1 = hh*k2[xmax-2]/12;
    double c2 = hh*k2[xmax-1]/12;
    for (int x=xmax-3; x>xc-3; x--) {
        y[x] = (y[x+1]*(2-10*c1) - y[x+2]*(1+c2)) / (1+c0);
        c2 = c1;
        c1 = c0;
        c0 = hh*k2[x-1]/12;
    }
}



inline double V_ho(double x)
{
    return x*x/2;
}



inline double V_lj(double x, double epsilon, double sigma)
{
    return 4*epsilon*(pow(sigma/x,12)-pow(sigma/x,6));
}


double E0(double h, double rmax, double (*V)(double))
{
    // we'll just assume that the minimum is near 0 ¯\_(ツ)_/¯
    return V(1e-2);
}


int nodeNumber(double * eigv, size_t N) // very rough
{
    int nodes = 0;
    int step = 10;
    for (int x=0; x<N; x+=step)
        if (eigv[x]*eigv[x-step] < 0) nodes++;
        
    return nodes;
}

void save2csv(Spectrum * spectra, int lmax, int nmax, int xmax)
{
    // save energy eigenvalues to csv file
    char filename[32];
    snprintf(filename, 32, "./eigenvalues.csv");
    FILE * eigenvalues;
    eigenvalues = fopen(filename, "w");
    if (eigenvalues == NULL)
        perror("error while writing on eigenvalues.csv");
    
    fprintf(eigenvalues, "l, n, E\n");
    for (int l=0; l<=lmax; l++)
        for (int n=0; n<nmax; n++)
            fprintf(eigenvalues, "%d, %d, %0.12f\n", l, n, spectra[l].EE[n]);
        
    // save eigenvectors to csv file
    FILE * eigenvectors;
    eigenvectors = fopen("./eigenvectors.csv", "w");
    if (eigenvectors == NULL)
        perror("error while writing on eigenvectors.csv");
    
    for (int l=0; l<=lmax; l++)
        for (int n=0; n<nmax; n++)
            fprintf(eigenvectors, "y%d%d,", l, n);
    fprintf(eigenvectors, "\n");
    for (int l=0; l<=lmax; l++)
        for (int n=0; n<nmax; n++)
            fprintf(eigenvectors, "%d,", l);
    fprintf(eigenvectors, "\n");
    for (int l=0; l<=lmax; l++)
        for (int n=0; n<nmax; n++)
            fprintf(eigenvectors, "%d,", n);
        
    fprintf(eigenvectors, "\n");
    
    for (int x=0; x<xmax; x++) {
        for (int l=0; l<=lmax; l++) {
            for (int n=0; n<nmax; n++)
                fprintf(eigenvectors, "%0.9f,", spectra[l].eigfuns[n*xmax+x]);
        }
        fprintf(eigenvectors, "\n");
    }
}



