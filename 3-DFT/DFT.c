#include "DFT.h"
//TODO 
// far funzionare numerov con vettore
// E_hartree
// calcolo energia cinetica TS
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
#if N==8
#define nmax 1
#else
#define nmax 2
#endif

#define rs rs_Na
#define rho_b  (3/(4*pi*rs*rs*rs))
#define Rc (pow(N,1/3)*rs)


int main(int argc, char** argv)
{ 
    // creates data folder and files to save data
    make_directory("Data");
    chdir("Data");
    
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
    double rho[gridlength];
    double rho_old[gridlength];
    for (int i=1; i < gridlength; i++)
    {
        rho[i] = pow(M_e, -(i*h)*(i*h))/(pi/2);
        rho_old[i] = rho[i];
    }
    double alpha = 0.1;
    double E1, E2; // total energy calculated in two different ways, to check for consistency
    Spectrum eig;
    eig.xmax = gridlength;
    eig.nmax = nmax;
    eig.EE = calloc(nmax, sizeof(double));
    eig.eigfuns = calloc(gridlength*nmax, sizeof(double));
    dArray bc; // boundary conditions
    bc.length = 2;
    bc.data = calloc(2, sizeof(double));
    bc.data[0] = 0;
    bc.data[1] = pow(rmax/gridlength,l);

    // check on the initial test density
    norm = normalizationFactor(rho, h, 0, gridlength);
    rho_ptr = &rho;
    printf("rho[0] = %f, rhoptr[0] = %f, rhonorm = %f\n", rho[0], rho_ptr[0], norm);
    
    
    
    /* Self-consistent iterative solution of the Kohn-Sham equation begins here */
    
    for (int t=0; t<MAXSTEPS; t++)
    {
        printf("\nITERATION %d\n", t);
        
        // solution of the kohn-sham equation with a trial density
        eig = numerov_var(nmax, 1, gridlength, rmax, h2m, 0.13, true, bc, rho, V_ks); //TODO rendere Vks un array o viceversa
        
        // calculation of the new density from the solutions
        for (int i=1; i < gridlength; i++)
        {
            rho[i] = ; // somma del modulo quadro delle autofunzioni trovate con numerov, per 2 di spin
            rho[i] =  alpha*rho[i] + (1-alpha)*rho_old[i]; // mixing of the solution with the old density (va messa prima o dopo?)
            rho_old[i] = rho_i;
        }

        // Consistency check through the energy
        E1 = E_ks(rho);
        E2 = sum(eig.EE, N) -0.5*+ E_xc(rho) // sum of the eigenvalues - 1/2 hartree energy - exchange 
        printf("E_ks = %f\tE_ = %f\tdiff = %f\n", E1, E2, E2-E1);
        
        // Check on the convergence by looking at how different is the new density
        delta = fabs(rho[0] - rho_old[0]);
        for (int i=1; i<gridlength; i++)
            if (fabs(rho[i]-rho_old[i]) > delta)
                delta = fabs(rho[i]-rho_old[i]);
            
        printf("delta = %f\n", delta);
        
    }
        
        
    free(bc.data);
    free(eig.eigfuns);
    
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
    return V_ks(r, rho_ptr)
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

inline double local_energy(double rho)
{
    return -0.75*pow(3*rho/pi, 1/3) - 0.44/(0.78 + 3/(4*pi*rho));
}

// Energy functional
double E_ks(double *rho)
{
    double E_k=0, E_ext=0, E_h=0;

    return E_k + E_ext(rho) + E_h(rho) + E_xc(rho);
}

// External energy
double E_ext(double *rho)
{
    double rhoVext[gridlength];
    for (int i=0; i < gridlength; i++)
        rhoVext = rho[i]*V_ext(i*h);

    return simpson_integral(rhoVext, gridlength, h);
}

// Hartree energy TODO check it
double E_h(double *rho)
{
    double E_H = 0;
    double integrand[gridlength];
    
    for (int i=0; i < gridlength; i++)
        integrand[i] = V_h(i*h, rho)*rho[i];
    
    for (int i=1; i < gridlength-1; i+=2)
        E_h += h*(rho[i-1]*integrand[i-1] + 4.*rho[i]*integrand[i] + rho[i+1]*integrand[i+1])/3.;

    return E_h/2;
}

// Exchange-correlation energy
double E_xc(double *rho)
{
    double E_xc = 0;
    for (int i=1; i < gridlength-1; i+=2)
        E_xc += h*(rho[i-1]*local_energy(rho[i-1]) + 4.*rho[i]*local_energy(rho[i]) + rho[i+1]*local_energy(rho[i+1]))/3.;

    return E_xc;
}
    

