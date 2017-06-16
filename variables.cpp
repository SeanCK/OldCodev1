#include   <stdio.h>
#include   <stdlib.h>
#include   <math.h>
#include   <iostream>
#include   <complex>


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


/* Functions  *********************************************************** */


double dot(double *r1, double *r2){ // definition of dot product

    double x;
    int i;
    for (i=0, x = 0.0; i < N_bath; i++)
        x += r1[i]*r2[i];
    return x;
}

double asyEps = 0.4; // new xxxxxxxxxxxxxxxxxxxxxxxx

double gam(double *R){

    double x;
    int i;

    for (i = 0, x = 0.0; i < N_bath; i++)
        x += c[i]*R[i];
    x += asyEps;    // asymmetric spin boson
    return -x;
}


double Hb(double *R, double *P){

    double x;
    int i;

    /* Bath Hamiltonian */

    for (i = 0, x = 0.0; i < N_bath; i++)
        x += P[i]*P[i] - mww[i]*R[i]*R[i];

    return x*0.5;
}

void dgamma(double *R){ //No returned value

    int i;

    for (i = 0; i < N_bath; i++)
        dgam[i] = -c[i];
}

void Fb(double *R){ //No returned value

    /* Pure Bath Force Field */

    int i;
    double x;

    for (i= 0; i < N_bath; i++)
        f[i] = mww[i]*R[i];
}


void F1(double *R){ //No returned value

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

void F2(double *R){ // No returned value

    int i;
    double g,h;

    /* 11 force field */

    g = gam(R);
    h = g/sqrt(ddd4 + g*g);

    for (i = 0; i< N_bath; i++)
        f[i] = mww[i]*R[i] + h*c[i];
}

double dE(double *R){

    double g;

    /* Energy difference between adibiatic surface (E1 - E0) */

    g = gam(R);
    g *= 4.0*g; //Is it g*g?
    return (sqrt(ddd + g));
}

double G(double *R){

    double x,g;

    g = gam(R);
    if (fabs(g/delta) < 1.0e-7)
        return (g/delta);
    x = (-delta + sqrt(ddd + 4*g*g))/(2*g);

    return x;
}

void dd(double*R){ //No returned value

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
    for (i = 0; i < N_bath; i++)
        dhat[i] /= abs_d;
    /* Note may need to change the sign of dhat here

    /* dhat now is not normalized, and equal d_{10}(R)
       -------------------------------------------------------------- */
}

void (*force[4])(double *); //No returned value


void integ_step(double *r, double *v, double dt, int Sa){ // No returned value

    /* ********* Velocity Verlet ************************** */

    double y;
    int i;

    y = 0.5*dt*dt;
    for (i = 0; i < N_bath; i++)
        r[i] += dt*v[i] + y*f[i];
    y = 0.5*dt;
    for (i = 0; i < N_bath; i++)
        v[i] += y*f[i];
    force[Sa](r);
    for (i = 0; i < N_bath; i++)
        v[i] += y*f[i];
}


void bath_para(double eta, double w_max){ //No returned value

    /* Parameters for bath (corresponding to an ohmic spectral density) */

    int i;
    double w_0;

    w_0 = (1 - exp(-w_max))/N_bath;

    for (i = 0; i < N_bath; i++){
        m[i] = 1.0;
        w[i] = -log( 1-(i+1)*w_0 );
        c[i] = sqrt(eta*w_0*m[i])*w[i];
    }
}


double U(double *r, double *v, int Sa, double t){

    double  dE0, phase,dt,x1,x2,x3,v1,v2,v3;
    int Nsteps, i;

    /* ******** Adiabatic Propagator *********************  */

    force[Sa](r);
    dt = timestep;

    if (t <= timestep){
        dt = t;
        Nsteps = 1;
    }
    else{
        Nsteps = t/dt +1;
        dt = t/Nsteps;
    }

    if ( (Sa == 0) || (Sa == 3)){
        for (i = 0; i < Nsteps; i++){
            integ_step(r , v,  dt, Sa);
        }
        return 0.0;
    }
    phase = dE(r)/2;
    for (i = 0; i < Nsteps; i++){
        integ_step(r , v,  dt, Sa);
        phase += dE(r);
    }
    phase -=dE(r)/2;
    phase*= dt;

    if (Sa == 1)
        phase *= -1.0;

    return phase;

    /* For some reason before I was returing  - phase  */

    /*   note dE = E1 - E0 > 0 , w_10 = dE , corresponding to the 1 0 pair
     of surfaces, equivalent to Sa = 2    */
}


double alpha,BETA1,BETA2, GAMMA1, GAMMA2,*b;



/* A load of  pointers initial densities and observables **********
    makes code changes easier   */

double (*phi)(double*, double*, int);

double (*dens_init[4])(double *, double *, int );

double (*obs[4])(double *, double *, int );

double (*obs1[4])(double *, double *, int );

double (* www[6][4][4])();



double sina, cosa, cosb1,cosb2, sinb1,sinb2,cosg1,cosg2,sing1,sing2, ppower,de;




//    Transition Matrices ___________________________________________________________________________-


/* Q1 */
double x;

double wwa0_00(){
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_01(){
    x = -sina;
    return x*0.5;
}

double wwa0_02(){
    x = -sina;
    return x*0.5;
}

double wwa0_03(){
    x = 1.0 - cosa;
    return x*0.5;
}

double wwa0_10(){
    x =  sina;
    return x*0.5;
}

double wwa0_11(){
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_12(){
    x = -1.0 + cosa;
    return x*0.5;
}

double wwa0_13(){
    x = -sina;
    return  x*0.5;
}

double wwa0_20(){
    x =  sina;
    return x*0.5;
}

double wwa0_21(){
    x = -1.0 + cosa ;
    return x*0.5;
}

double wwa0_22(){
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_23(){
    x = -sina ;
    return x*0.5;
}

double wwa0_30(){
    x = 1.0 - cosa;
    return x*0.5;
}

double wwa0_31(){
    x =  sina;
    return x*0.5;
}

double wwa0_32(){
    x =  sina;
    return x*0.5;
}

double wwa0_33(){
    x = 1.0 + cosa;
    return x*0.5;
}

/* _____________________________________________  */

//         W_{a 1}

double wwa1_00(){
    return 9999.0;
}


double wwa1_01(){
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_02(){
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_03(){
    x = Pdotdhat*Pdotdhat - 2.0*de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_10(){
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_11(){
    return 9999.0;
}

double wwa1_12(){
    return 9999.0;
}

double wwa1_13(){
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_20(){
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_21(){
    return 9999.0;
}

double wwa1_22(){
    return 9999.0;
}

double wwa1_23(){
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_30(){
    x = Pdotdhat*Pdotdhat + 2.0*de;
    return sqrt(x);
}

double wwa1_31(){
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_32(){
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_33(){
    return 9999.0;
}


//   W_{a2}

/* _____________________________________________  */

double wwa2_00(){
    x = sina -  sinb2;
    return x;
}

double wwa2_01(){
    x = cosb2;
    return x;
}

double wwa2_02(){
    x = cosb2;
    return x;
}

double wwa2_03(){
    x = sina  + sinb2;
    return  x;
}

double wwa2_10(){
    x = cosa;
    return x;
}

double wwa2_11(){
    return 0.0;
}

double wwa2_12(){
    return 0.0;
}

double wwa2_13(){
    x = cosa;
    return x;
}

double wwa2_20(){
    x = cosa;
    return x;
}

double wwa2_21(){
    return 0.0;
}

double wwa2_22(){
    return 0.0;
}

double wwa2_23(){
    x = cosa;
    return x;
}

double wwa2_30(){
    x = sina +  sinb2;
    return  -x;
}

double wwa2_31(){
    x = cosb2;
    return x;
}

double wwa2_32(){
    x = cosb2;
    return x;
}

double wwa2_33(){
    x = -sina  + sinb2;
    return  x;
}


//   Wb0

double wwb0_00(){
    x = cosa;
    return x*x;
}

double wwb0_01(){
    x = cosa*sina;
    return x;
}

double wwb0_02(){
    x =  cosa*sina;
    return x;
}

double wwb0_03(){
    x = sina;
    return 1.0 + x*x;
}

double wwb0_10(){
    x = -cosa*sina;
    return x;
}

double wwb0_11(){
    x = cosa;
    return x*x;
}

double wwb0_12(){
    x = cosa;
    return x*x;
}

double wwb0_13(){
    x = cosa*sina;
    return -x;
}

double wwb0_20(){
    x = cosa*sina;
    return -x;
}

double wwb0_21(){
    x = cosa;
    return x*x;
}

double wwb0_22(){
    x = cosa;
    return x*x;
}

double wwb0_23(){
    x = cosa*sina;
    return  x;
}

double wwb0_30(){
    x = sina;
    return 1.0 + x*x;
}

double wwb0_31(){
    x = cosa*sina;
    return -x;
}

double wwb0_32(){
    x = cosa*sina;
    return -x;
}

double wwb0_33(){
    x = cosa;
    return  x*x;
}



//  W_{b1}


double wwb1_00(){
    x = 1.0 - sina*sing1;
    return x;
}

double wwb1_01(){
    x = cosg1*sina;
    return x;
}

double wwb1_02(){
    x = cosg1*sina;
    return x;
}

double wwb1_03(){
    x = 1.0 + sina*sing1;
    return x;
}

double wwb1_10(){
    x = cosa*sing1;
    return -x;
}

double wwb1_11(){
    x = cosa*cosg1;
    return x;
}

double wwb1_12(){
    x = cosa*cosg1;
    return x;
}

double wwb1_13(){
    x = cosa*sing1;
    return x;
}

double wwb1_20(){
    x = cosa*sing1;
    return -x;
}

double wwb1_21(){
    x = cosa*cosg1;
    return x;
}

double wwb1_22(){
    x = cosa*cosg1;
    return x;
}

double wwb1_23(){
    x = cosa*sing1;
    return x;
}

double wwb1_30(){
    x = 1.0 + sina*sing1;
    return x;
}

double wwb1_31(){
    x = -cosg1*sina;
    return x;
}

double wwb1_32(){
    x = -cosg1*sina;
    return x;
}

double wwb1_33(){
    x = 1.0 - sina*sing1;
    return x;
}


//  W_{b2}


double wwb2_00(){
    x = 1.0 - sina*sing2;
    return x;
}


double wwb2_01(){
    x = cosg2*sina;
    return x;
}

double wwb2_02(){
    x = cosg2*sina;
    return x;
}

double wwb2_03(){
    x = 1.0 + sina*sing2;
    return x;
}

double wwb2_10(){
    x = cosa*sing2;
    return -x;
}

double wwb2_11(){
    x = cosa*cosg2;
    return x;
}

double wwb2_12(){
    x = cosa*cosg2;
    return x;
}

double wwb2_13(){
    x = cosa*sing2;
    return x;
}

double wwb2_20(){
    x = cosa*sing2;
    return -x;
}

double wwb2_21(){
    x = cosa*cosg2;
    return x;
}

double wwb2_22(){
    x = cosa*cosg2;
    return x;
}

double wwb2_23(){
    x = cosa*sing2;
    return x;
}

double wwb2_30(){
    x = 1.0 + sina*sing2;
    return x;
}

double wwb2_31(){
    x = cosg2*sina;
    return -x;
}

double wwb2_32(){
    x = cosg2*sina;
    return -x;
}

double wwb2_33(){
    x = 1.0 - sina*sing2;
    return x;
}


void setwww(){

    // W_{a0}

    www[0][0][0] = wwa0_00;
    www[0][0][1] = wwa0_01;
    www[0][0][2] = wwa0_02;
    www[0][0][3] = wwa0_03;
    www[0][1][0] = wwa0_10;
    www[0][1][1] = wwa0_11;
    www[0][1][2] = wwa0_12;
    www[0][1][3] = wwa0_13;
    www[0][2][0] = wwa0_20;
    www[0][2][1] = wwa0_21;
    www[0][2][2] = wwa0_22;
    www[0][2][3] = wwa0_23;
    www[0][3][0] = wwa0_30;
    www[0][3][1] = wwa0_31;
    www[0][3][2] = wwa0_32;
    www[0][3][3] = wwa0_33;

    // W_{a1}

    www[1][0][0] = wwa1_00;
    www[1][0][1] = wwa1_01;
    www[1][0][2] = wwa1_02;
    www[1][0][3] = wwa1_03;
    www[1][1][0] = wwa1_10;
    www[1][1][1] = wwa1_11;
    www[1][1][2] = wwa1_12;
    www[1][1][3] = wwa1_13;
    www[1][2][0] = wwa1_20;
    www[1][2][1] = wwa1_21;
    www[1][2][2] = wwa1_22;
    www[1][2][3] = wwa1_23;
    www[1][3][0] = wwa1_30;
    www[1][3][1] = wwa1_31;
    www[1][3][2] = wwa1_32;
    www[1][3][3] = wwa1_33;

    // W_{a2}

    www[2][0][0] = wwa2_00;
    www[2][0][1] = wwa2_01;
    www[2][0][2] = wwa2_02;
    www[2][0][3] = wwa2_03;
    www[2][1][0] = wwa2_10;
    www[2][1][1] = wwa2_11;
    www[2][1][2] = wwa2_12;
    www[2][1][3] = wwa2_13;
    www[2][2][0] = wwa2_20;
    www[2][2][1] = wwa2_21;
    www[2][2][2] = wwa2_22;
    www[2][2][3] = wwa2_23;
    www[2][3][0] = wwa2_30;
    www[2][3][1] = wwa2_31;
    www[2][3][2] = wwa2_32;
    www[2][3][3] = wwa2_33;

    //  wb0

    www[3][0][0] = wwb0_00;
    www[3][0][1] = wwb0_01;
    www[3][0][2] = wwb0_02;
    www[3][0][3] = wwb0_03;
    www[3][1][0] = wwb0_10;
    www[3][1][1] = wwb0_11;
    www[3][1][2] = wwb0_12;
    www[3][1][3] = wwb0_13;
    www[3][2][0] = wwb0_20;
    www[3][2][1] = wwb0_21;
    www[3][2][2] = wwb0_22;
    www[3][2][3] = wwb0_23;
    www[3][3][0] = wwb0_30;
    www[3][3][1] = wwb0_31;
    www[3][3][2] = wwb0_32;
    www[3][3][3] = wwb0_33;

    // W_{b1}

    www[4][0][0] = wwb1_00;
    www[4][0][1] = wwb1_01;
    www[4][0][2] = wwb1_02;
    www[4][0][3] = wwb1_03;
    www[4][1][0] = wwb1_10;
    www[4][1][1] = wwb1_11;
    www[4][1][2] = wwb1_12;
    www[4][1][3] = wwb1_13;
    www[4][2][0] = wwb1_20;
    www[4][2][1] = wwb1_21;
    www[4][2][2] = wwb1_22;
    www[4][2][3] = wwb1_23;
    www[4][3][0] = wwb1_30;
    www[4][3][1] = wwb1_31;
    www[4][3][2] = wwb1_32;
    www[4][3][3] = wwb1_33;

    // W_{b2}

    www[5][0][0] = wwb2_00;
    www[5][0][1] = wwb2_01;
    www[5][0][2] = wwb2_02;
    www[5][0][3] = wwb2_03;
    www[5][1][0] = wwb2_10;
    www[5][1][1] = wwb2_11;
    www[5][1][2] = wwb2_12;
    www[5][1][3] = wwb2_13;
    www[5][2][0] = wwb2_20;
    www[5][2][1] = wwb2_21;
    www[5][2][2] = wwb2_22;
    www[5][2][3] = wwb2_23;
    www[5][3][0] = wwb2_30;
    www[5][3][1] = wwb2_31;
    www[5][3][2] = wwb2_32;
    www[5][3][3] = wwb2_33;

}



/* ************** Observables and initial density Matrices *********  */


/* Definition of initial density matrix element */
double z;
double *y;
double g = G(y);
double gg = g*g;

double wigner_harm_osc(double *y, double *p){ /* set wigner to one has it drops out of calculation for our initial
    density */

    double prod,y1,y2;
    int i;

    return 1.0;
}


double dens_init_0(double *y,double *p, int k){
    z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);
    return (z*wigner_harm_osc(y,p));
}

double dens_init_1(double *y,double *p, int k){ //=2
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(y,p));
}

double dens_init_2(double *y,double *p, int k){ //=1
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(y,p));
}

double dens_init_3(double *y,double *p, int k){
    z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);
    return (z*wigner_harm_osc(y,p));
}


double obs_0(double *y,double *p, int k){ //=-3
    z = 2.0*g/(1 + gg);
    return z;
}

double obs_1(double *y,double *p, int k){ //=2
    z = (gg-1)/(1 + gg);
    return z;
}

double obs_2(double *y,double *p, int k){ //=1
    z =  (gg-1)/(1 + gg);
    return z;
}

double obs_3(double *y,double *p, int k){ //=0
    z = -2.0*g/(1 + gg);
    return z;
}


/* These are the matrix elements of the Hamiltonian */

double H_0(double *y,double *p, int k){
    z = Hb(y,p) - dE(y)*0.5;
    return z;
}

double H_1(double *y,double *p, int k){ //=2
    return 0.0;
}

double H_2(double *y,double *p, int k){ //=1
    return 0.0;
}

double H_3(double *y,double *p, int k){
    z = Hb(y,p) + dE(y)*0.5;
    return z;
}
