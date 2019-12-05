#include "DFT.h"
//TODO 
// calcolo energia cinetica TS (serve a qualcosa?)
// vari printf 
// usare save2csv
// capire boundaries
// capire l

#define pi M_PI
#define M_e 2.718281828459045
#define h2m 0.5
#define rmax 5.0
#define gridlength 5000
#define h (rmax/gridlength)
#define MAXSTEPS 1000

#define rs_Na 3.93
#define rs_K 4.86
#define rho_bNa 0.00393309
#define rho_bK 0.00207971

#define N 8
#define NMAX 1 // should be 1 for Z=8 and 2 for Z=20???
#define LMAX 1 // 1 (p) for Z=8 and 2 (d) for Z=20???

#define rs rs_Na
#define rho_b  (3/(4*pi*rs*rs*rs))
#define Rc (pow(N,1/3)*rs)

double * RHO;

int main(int argc, char** argv)
{ 
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
    FILE * ksenergy; // saves energy from the two different calculations
    ksenergy = fopen("ksenergy.csv", "w");
    FILE * ksdensity; // saves densities, mixed and unmixed
    ksdensity = fopen("ksdensity.csv", "w");

    
    /*int N;
    if (argc == 2)
    {
        N = (int)strtol(argv[1], NULL, 10);
    }
    else    {
        // asks user for Z
        printf("Enter Z (8 or 20): ");
        scanf("%d",&N);
    }
    if (N==8) 
        double rs = rs_Na;
    else 
        double rs = rs_K;
    
    static double rho_b = 3/(4*pi*rs*rs*rs);
    static double Rc = pow(N,1/3)*rs;*/
    
    double norm, delta;
    //double rho = calloc(gridlength, sizeof(double));
    RHO = calloc(gridlength, sizeof(double));
    double * rho_old = calloc(gridlength, sizeof(double));
    // starting guess for the density is a normalized gaussian
    for (int i=1; i < gridlength; i++)
    {
        RHO[i] = pow(M_e, -(i*h)*(i*h))/(pi/2);
        rho_old[i] = RHO[i];
    }
    double alpha = 0.1; // mixing coefficient of the densities
    double E1, E2; // total energy calculated in two different ways, to check for consistency
    
    Spectrum spectra[LMAX+1];
    for (int l=0; l<=LMAX; l++) {
        spectra[l].xmax = gridlength;
        spectra[l].nmax = NMAX;
        spectra[l].EE = calloc(NMAX, sizeof(double));
        spectra[l].eigfuns = calloc(gridlength*NMAX, sizeof(double));
    }
    dArray bc; // boundary conditions in 0
    bc.length = 2;
    bc.data = calloc(2, sizeof(double));
    bc.data[0] = 0.0; //NOTE should check this

    // check on the initial test density
    norm = normalizationFactor(RHO, h, 0, gridlength);
    printf("rho[0] = %f, rhonorm = %f\n", RHO[0], norm);
    
    
    
    /* Self-consistent iterative solution of the Kohn-Sham equation begins here */
    
    for (int t=0; t<MAXSTEPS; t++)
    {
        printf("\nITERATION %d\n", t);
        
        // solution (Kohn-Sham orbitals) of the kohn-sham equation with a trial density RHO, for different angular momentum values
        for (int l=0; l<=LMAX; l++) {
            bc.data[1] = pow(rmax/gridlength,l);
            spectra[l] = numerov(NMAX, 1, gridlength, rmax, h2m, 0.13, true, bc, V_ks_rho);
        }
        save2csv(spectra, LMAX, NMAX, gridlength);
        
        // calculation of the new density from the solutions, including the mix with the old density
        for (int i=1; i < gridlength; i++)
        {
            // somma del modulo quadro delle autofunzioni trovate con numerov, times degeneration
            if (N==8)   {
                // adds the 2 electrons in the 1s orbital + the 6 in the 1p
                RHO[i] = 2*spectra[0].eigfuns[i]*spectra[0].eigfuns[i] + 6*spectra[1].eigfuns[i]*spectra[1].eigfuns[i];
            }
            else if (N==20)  {
                // adds the 2 electrons in the 1s orbital + the 6 in the 1p
                RHO[i] = 2*spectra[0].eigfuns[i]*spectra[0].eigfuns[i] + 6*spectra[1].eigfuns[i]*spectra[1].eigfuns[i];
            }
            fprintf(ksdensity, ",y%d%d", l, n);
            RHO[i] =  alpha*RHO[i] + (1-alpha)*rho_old[i]; // mixing of the solution with the old density (va messa prima o dopo ks?)
            rho_old[i] = RHO[i];
        }

        // Consistency check through the energy
        E1 = E_ks(RHO);
        E2 = spectra[0].EE[0]*2 + spectra[1].EE[0]*6 - E_H(RHO) + E_XC(RHO); // sum of the eigenvalues - 1/2 hartree energy - exchange 
        printf("E_ks = %f\tE_ = %f\tdiff = %f\n", E1, E2, E2-E1);
        
        // Check on the convergence by looking at how different is the new density
        delta = fabs(RHO[0] - rho_old[0]);
        for (int i=1; i<gridlength; i++)
            if (fabs(RHO[i]-rho_old[i]) > delta)
                delta = fabs(RHO[i]-rho_old[i]);
            
        printf("delta = %f\n", delta);
        
    }
        
        
    free(bc.data);
    for (int l=0; l<=LMAX; l++)
        free(spectra[l].eigfuns);
    
    return 0;
}



// Kohn-Sham potential
double V_ks(double r, const double *rho)
{    
    // Exchange-correlation potential
    double V_xc = (-0.25*pow(3/pi,1./3)*pow(rho[(int)(h*r)],-2./3) - 0.14667*pow(4/(3*pi),1./3)*pow(rho[(int)(h*r)],-2./3)/(0.78+3/(4*pi*rho[(int)(h*r)])))*local_energy(rho[(int)(h*r)]) + local_energy(rho[(int)(h*r)]);
    
    return V_ext(r) + V_h(r, rho) + V_xc;
}

double V_ks_rho(double r)
{
    return V_ks(r, RHO);
}

double V_ext(double r)
{
    if (r>Rc)
        return -4*pi*rho_b*Rc*Rc*Rc/(3*r);
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc);
}

// Hartee potential in the Kohn-Sham equation, also used to compute the Hartree energy
double V_h(double r, const double *rho)
{
    double V_h = 0.;
    for (int i=1; i < (int)(h*r)-1; i+=2)
        V_h += h * (rho[i-1]*(h*(i-1)) + 4.*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3.;
    
    for (int i=(int)(h*r); i < gridlength-1; i+=2)
        V_h += h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4.*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / (3.*r);
    
    return 4*pi*V_h;
}



// Energy functional
double E_ks(double *rho)
{
    return T_S(rho) + E_ext(rho) + E_H(rho) + E_XC(rho);
}

// External energy
double E_ext(double *rho)
{
    double rhoVext[gridlength];
    for (int i=0; i < gridlength; i++)
        rhoVext[i] = rho[i]*V_ext(i*h);

    return simpson_integral(rhoVext, gridlength, h);
}

// Hartree energy TODO check it bc i'm really not sure why i'm programming at BUC 10 minutes before closure
double E_H(double *rho)
{
    double E_h = 0;
    double integrand[gridlength];
    
    for (int i=0; i < gridlength; i++)
        integrand[i] = V_h(i*h, rho)*rho[i];
    
    for (int i=1; i < gridlength-1; i+=2)
        E_h += h*(rho[i-1]*integrand[i-1] + 4.*rho[i]*integrand[i] + rho[i+1]*integrand[i+1])/3.;

    return E_h/2; //check this
}

// Exchange-correlation energy
double E_XC(double *rho)
{
    double E_xc = 0;
    for (int i=1; i < gridlength-1; i+=2)
        E_xc += h*(rho[i-1]*local_energy(rho[i-1]) + 4.*rho[i]*local_energy(rho[i]) + rho[i+1]*local_energy(rho[i+1]))/3.;

    return E_xc;
}

// Kohn–Sham kinetic energy (should take the ks orbitals as input)
double T_S(double *phi)
{
    double integrand[gridlength-4]; //excludes extrema used as boundary conditions
    
    for(int i=2; i<gridlength-2; i++)
        integrand[i-2] = der5(phi,i,h) * der5(phi,i,h);
    
    return h2m * simpson_integral(integrand, gridlength-4, h);
}


// Local energy used in the LDA
inline double local_energy(double rho)
{
    return -0.75*pow(3*rho/pi, 1/3) - 0.44/(0.78 + 3/(4*pi*rho));
}
