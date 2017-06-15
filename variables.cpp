#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   <numeric>


using namespace std;

#include <gsl/gsl_rng.h>

#define PI 3.141592653589793


/* Declaration of external/global variables and pointers    ****************** */

char  datafilename[80];
FILE   *stream;
const gsl_rng_type * TT;
gsl_rng * rr;

int rng_seed;
int N_bath,MAX_tr,tseed;
int *SS, *jump, rand_m, LL;
int rho_top_init[4], obs_top[4];
int  N_slice,Nsample,Ncut,**hist;
int  rngstreamnum, nrngstreams, *rngstream;
double *m,*c,*w,*d,delta,beta, **RR, **PP,  *dhat, *f, *dgam,*Pperp;
double *mww;
double *mu,*sig,ddd4,ddd,*dtdtm,ranVector[10001];
double *meann;
double abs_d, norm_c, timestep, sigma, xinit, loco1,TSLICE,Dt,yoyo = .2,Pdotdhat, nojumpknob = 1.0E-6;
double  *abszsum0, *abszsum1, *argzsum0, *argzsum1, **realsum, **imagsum;
double  *habszsum0, *habszsum1, *hargzsum0, *hargzsum1, **hrealsum, **himagsum;
complex<double> I(0,1);

int t_strobe, Nblock = 1024; /* t_strobe is the frequency at which results for slices are printed,
                                 Nblock is the size of sub-ensembles */




/* definition of dot product  *********************************************************** */


double dot(double *r1, double *r2){

    double x;
    int i;

    for (i=0, x = 0.0; i < N_bath; i++){
        x += r1[i]*r2[i];
	}
    return x;

}

double asyEps = 0.4; // new xxxxxxxxxxxxxxxxxxxxxxxx

double gam(double *R){

    double x;
    int i;

    for (i = 0, x = 0.0; i < N_bath; i++){
        x += c[i]*R[i];
	}
    x += asyEps;    // asymmetric spin boson

    return -x;
}


double Hb(double *R, double *P){

    double x;
    int i;

    /* Bath Hamiltonian */

    for (i = 0, x = 0.0; i < N_bath; i++){
        x += P[i]*P[i] - mww[i]*R[i]*R[i];
	}
    return x*0.5;
}


void dgamma(double *R){

    int i;

    for (i = 0; i < N_bath; i++){
        dgam[i] = -c[i];
	}
}

void Fb(double *R){

    /* Pure Bath Force Field */

    int i;
    double x;

    for (i= 0; i < N_bath; i++){
        f[i] = mww[i]*R[i];
	}
}


void F1(double *R){

    int i;
    double g,h;

    /* 00 force field   */

    g = gam(R);
    h = g/sqrt(ddd4 + g*g);
    for (i = 0; i< N_bath; i++){
        f[i]  = mww[i]*R[i] -  h*c[i];
    }
    /* ok  note mww[i] = - m*w[i]*w[i] as calculated in main */
}


void F2(double *R){

    int i;
    double g,h;

    /* 11 force field */

    g = gam(R);
    h = g/sqrt(ddd4 + g*g);

    for (i = 0; i< N_bath; i++){
        f[i] = mww[i]*R[i] + h*c[i];
	}
}


double dE(double *R){

    double g;

    /* Energy difference between adibiatic surface (E1 - E0) */


    g = gam(R);
    g *= 4.0*g;
    return (sqrt(ddd + g));

    /* This is E1 - E0   ok */
}


double G(double *R){

    double x,g;

    g = gam(R);
    if (fabs(g/delta) < 1.0e-7)
        return (g/delta);
    x = (-delta + sqrt(ddd + 4*g*g))/(2*g);

    return x;

}


void dd(double*R){

    double x1,x2,x3;
    int i;

    /* Computing d = <1|d/dr|0>  */

    x2 = gam(R);
    dgamma(R);
    if ( fabs(x2) < 1.0e-4)
        x3 = 1/delta;
    else {
        x1 = G(R);
        x3 = -x1/x2 + 2.0/(delta + 2.0*x2*x1);
        x3 = x3/(1.0 + x1*x1);
    }
    for (i = 0,abs_d = 0.0; i < N_bath; i++){
        dhat[i] = -dgam[i]*x3;   // 1-3-05 put - sign in
        abs_d += dhat[i]*dhat[i];
    }
    abs_d = sqrt(abs_d);
    for (i = 0; i < N_bath; i++) {
        dhat[i] /= abs_d;
        /* Note may need to change the sign of dhat here*/
        }
        /* dhat now is not normalized, and equal d_{10}(R)
           -------------------------------------------------------------- */
    }



/* ************** Observables and initial density Matrices *********  */

/* A load of pointers initial densities and observables ********** 
	makes code changes easier   */
double (*phi)(double*, double*, int);
double (*dens_init[4])(double*, double*, int);
double (*obs[4])(double*, double*, int);
double (*obs1[4])(double*, double*, int);
double (*www[6][4][4])();



/* Definition of initial density matrix element */


double wigner_harm_osc(double *x, double *p){
    /* set wigner to one has it drops out of calculation for our initial
    density */
    return 1.0;
}


/*Initial Densities */
/*Can all be parallized*/


double z;
double g=G(x);
double gg=g*g;


double dens_init_0(double *x,double *p, int k){
	z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);
	return (z*wigner_harm_osc(x,p));
}


double dens_init_1(double *x,double *p, int k){
	z = 0.5*(gg - 1.0)/(1 + gg);
	return (z*wigner_harm_osc(x,p));
}

double dens_init_2(double *x,double *p, int k){
	z = 0.5*(gg - 1.0)/(1 + gg);
	return (z*wigner_harm_osc(x,p));
}

double dens_init_3(double *x,double *p, int k){
	z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);
	return (z*wigner_harm_osc(x,p));
}



double obs_0(double *x,double *p, int k){
	z = 2.0*g/(1 + gg);
	return z;
}

double obs_1(double *x,double *p, int k){
	z = (gg-1)/(1 + gg);
	return z;
}

double obs_2(double *x,double *p, int k){
	z =  (gg-1)/(1 + gg);
	return z;
}

double obs_3(double *x,double *p, int k){
    z = -2.0*g/(1 + gg);
    return z;
}



/* These are the matrix elements of the Hamiltonian */

double H_0(double *x,double *p, int k){
	z = Hb(x,p) - dE(x)/2.0;
	return z;
}

double H_1(double *x,double *p, int k){
	return 0.0;
}

double H_2(double *x,double *p, int k){
	return 0.0;
}

double H_3(double *x,double *p, int k){
	z = Hb(x,p) + dE(x)/2.0;
	return z;
}
