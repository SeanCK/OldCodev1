#ifndef VARIABLE
#define VARIABLE

#include <gsl/gsl_rng.h>


int N_bath;
int  N_slice, Nsample, Ncut;

const gsl_rng_type * TT;
gsl_rng * rr;

double (*dens_init[4])(double *, double *, int );
double (*obs[4])(double *, double *, int );
double (*obs1[4])(double *, double *, int );
double (* www[6][4][4])();
void (*force[4])(double *);


int *SS, **hist;
double *mww;
double *meann, *sig;
double ddd4, ddd, ppower, de;
double *m, *c, *w, *d, delta, **RR, **PP, *dhat, *f, *dgam, *Pperp;
double abs_d, timestep, TSLICE, Dt, Pdotdhat;

double sina, cosa, alpha;

double  *abszsum0, *abszsum1, *argzsum0, *argzsum1, **realsum, **imagsum;
double  *habszsum0, *habszsum1, *hargzsum0, *hargzsum1, **hrealsum, **himagsum;

#endif