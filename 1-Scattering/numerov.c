// TODO
/*
 * cosa serve cambiare per scattering
 * migliorare find_zero in modo che pigli 0 appena dopo xc?
 * oppure passare a Numerov funzione precalcolata per tipo di potenziale?
 * aggiungere stocasticità a mini gradient-descent in E0
*/

#include "numerov.h"
#include "matematicose.c"
#include "misccose.c"

// hbar/2m, to change to LJ natural units at point 4
#define h2m 0.5


// Performs the whole algorithm and finds the spectrum up to the n-th level
// Returns a Spectrum struct
Spectrum numerov(int nmax, int l, int xmax, double rmax, double Estep, bool normalize, double (*V)(double))
{
    Spectrum sp;
    sp.EE = calloc(nmax, sizeof(double));
    sp.eigfuns = calloc(xmax*nmax, sizeof(double));
    double h = rmax/xmax;
    double xc_float; int xc; // point corresponding to the classical barrier
    int nfound = 0; // number of eigenvalues found
    double E = E0(V,h,rmax)+V(1e-1); // passare come argomento a Numerov?
    double E1, E2, E_;
    double delta, delta1, delta2, prevdelta=0;  // difference in the forward and backward log-derivatives
    double prev_yc = 0; // needed because when function flips sign the delta also changes sign producing false positives
    double V_[xmax], centrifugal[xmax], k2[xmax];
    double yf[xmax], yb[xmax];  // calculated functions, respectively before xc and after xc
    
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
        xc_float = findzero_last(V, E, 0., rmax, -h/4, h/4); // finds the 0 of V(x)-E with the secant method.
        //printf("xc_float %f, vs ho %f\n",xc_float,sqrt(2*E));
        //xc_float = sqrt(2*E); // !!!! vale solo per armonico !!!!
        xc = (int)round(xc_float/h);
        
        for (int x=0; x<xmax; x++)  
            k2[x] = (E-V_[x])/h2m - centrifugal[x];
        
        // Boundary conditions at rmax
        yb[xmax-1] = exp(-sqrt(fabs(E)/h2m)*xmax*h);
        yb[xmax-2] = exp(-sqrt(fabs(E)/h2m)*(xmax-1)*h);
        
        numerov_forward(h*h, xc, k2, yf);
        numerov_backward(h*h, xc, xmax, k2, yb);
        delta = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc];
        //printf("yf[xc-1] = %f,  yf[xc] = %f,  yf[xc+1]=%f\n", yf[xc-1], yf[xc], yf[xc+1]);
        //printf("yb[xc-1] = %f,  yb[xc] = %f,  yb[xc+1]=%f\n", yb[xc-1], yb[xc], yb[xc+1]);
        //printf("derforward = %f,  derback = %f\n", der3(yf,xc,h), der3(yb,xc,h)/yb[xc]);
        //printf("E = %f\tDelta rough = %f\n", E, delta);
        
        // If there is a change in sign, start the finer search of the 0 of delta(E) with the secant method
        if (delta*prevdelta < 0 && yf[xc]*prev_yc > 0)
        {
            printf("\nFound a point of inversion at %f - %f\n", E, E-Estep);
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
                if (E2 > V_[xmax-1]) perror("\n!!! Energy has exploded over limits !!! \n");
                if (E2 < 0) perror("\n!!! Energy is negative !!! \n");
                
                // Calculation of the next delta2 value
                xc_float = findzero_last(V, E2, 0., rmax, -h/4, h/4); // finds the 0 of V(x)-E with the secant method.
                //xc_float = sqrt(2*E2); // !!!! vale solo per armonico !!!!
                xc = (int)round(xc_float/h);
                
                for (int x=0; x<xmax; x++)
                    k2[x] = (E2-V_[x])/h2m - centrifugal[x];
                
                yb[xmax-1] = exp(-sqrt(fabs(E)/h2m)*xmax*h);
                yb[xmax-2] = exp(-sqrt(fabs(E)/h2m)*(xmax-1)*h);
                numerov_forward(h*h, xc, k2, yf);
                numerov_backward(h*h, xc, xmax, k2, yb);
                delta2 = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc];
            }
            
            printf("E%d = %.9f\n\n", nfound, E2);
            // Saves the eigenstate found
            sp.EE[nfound] = E2;
            for (int x=0; x<xc; x++)
                sp.eigfuns[xmax*nfound + x] = yf[x];
            for (int x=xc; x<xmax; x++)
                sp.eigfuns[xmax*nfound + x] = yb[x] * yf[xc]/yb[xc]; //impose continuity
            
            if (normalize==true)    {
                double norm = normalizationFactor(sp.eigfuns, h, xmax*nfound, xmax*(nfound+1));
                for (int x = xmax*nfound; x < xmax*(nfound+1); x++)
                    sp.eigfuns[x] = sp.eigfuns[x] / norm;
            }
                
            nfound++;
        }
        
        prevdelta = delta;
        prev_yc = yf[xc];
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

inline double V_fastlj(double dr2)
{
     double dr6 = dr2*dr2*dr2;
     return 4*(1.0/(dr6*dr6) - 1.0/dr6);
}


double E0_stupid(double (*V)(double), double h, double rmax)
{
    // we'll just assume that the minimum is near 0 ¯\_(ツ)_/¯
    return V(1e-1);
}

// Very basic gradient descent toward the potential minimum (currently the nearest minimum to the middle of the grid)
double E0(double (*V)(double), double h, double rmax)
{
    double r = rmax/2; // starting point
    double scale = fabs(V(rmax)-V(r)); // to get an order of magnitude of the variation of V 
    double gamma = scale/500;
    double grad = 10.;
    
    while (fabs(grad) > 1e-5)
    {
        grad = der5_c(V,r,h);
        r -= gamma*grad;
        //printf("grad = %f\tr=%f\n", grad, r);
    }
    
    return V(r);
}


int nodeNumber(const double * eigv, size_t N) // very rough
{
    int nodes = 0;
    int step = 10;
    for (int x=0; x<N; x+=step)
        if (eigv[x]*eigv[x-step] < 0) nodes++;
        
    return nodes;
}

double normalizationFactor(const double * eigv, double h, int x1, int x2)
{
    double sqmod[x2-x1];
    
    for(int i=x1; i<x2; i++)
        sqmod[i-x1] = eigv[i]*eigv[i];
    
    return integrale_simpson(sqmod, x2-x1, h);
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
    
    fprintf(eigenvectors, "y%d%d", 0, 0);
    for (int l=0; l<=lmax; l++) {
        for (int n=0; n<nmax; n++)  {
            if (!(l==0 && n==0))
                fprintf(eigenvectors, ",y%d%d", l, n);
        }
    }
    fprintf(eigenvectors, "\n0");
    for (int l=0; l<=lmax; l++) {
        for (int n=0; n<nmax; n++)  {
            if (!(l==0 && n==0))
                fprintf(eigenvectors, ",%d", l);
        }
    }
    fprintf(eigenvectors, "\n0");
    for (int l=0; l<=lmax; l++) {
        for (int n=0; n<nmax; n++)  {
            if (!(l==0 && n==0))
                fprintf(eigenvectors, ",%d", n);
        }
    }
        
    fprintf(eigenvectors, "\n");
    
    for (int x=0; x<xmax; x++) {
        fprintf(eigenvectors, "%0.9f", spectra[0].eigfuns[0*xmax+x]);
        for (int l=0; l<=lmax; l++) {
            for (int n=0; n<nmax; n++)  {
                if (!(l==0 && n==0))
                    fprintf(eigenvectors, ",%0.9f", spectra[l].eigfuns[n*xmax+x]);
            }
        }
        fprintf(eigenvectors, "\n");
    }
}



