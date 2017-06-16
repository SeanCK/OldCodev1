#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "dmatrix.c"
#include   "imatrix.c"
#include   "omp.h"
using namespace std;
/* #define SIMPLE_SPRNG		 simple interface                         */
// #include "sprng.h"
/* SPRNG header file     see NSCA random number generator   library   *** */
#include <gsl/gsl_rng.h>
/*
#define SEED 985456376
*/

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

    for (i=0, x = 0.0; i < N_bath; i++)
        x += r1[i]*r2[i];

    return x;

}

/* Random number stuff ***************************************************************  */


double  gauss1(double mean_x, double sigma_x,int i){

    double x1,y1,y2,y3;

    /* The tree here is to avoid a minor instability ***************************       */

    y1 = ranVector[i];

    if (fabs(y1) < 1.0e-200){
        y1 = gsl_rng_uniform (rr);
        if (fabs(y1) < 1.0e-200)
            y1 = (gsl_rng_uniform (rr));
        if (fabs(y1) < 1.0e-200)
            y1 = (gsl_rng_uniform (rr));
        if (fabs(y1) < 1.0e-200)
            y1 = gsl_rng_uniform (rr);
    }

    y2 = ranVector[i+N_bath];
    y3 = sqrt(-2*log(y1));
    x1 = y3*cos(2*PI*y2);
    return (mean_x =sigma_x*x1);

}

double dE(double *R);

void randnums(int rand_dim, double  *rand_vec){

    int i;

    for (i = 0; i < rand_dim; i++){
        rand_vec[i] = gsl_rng_uniform (rr);

    }
}


void gauss_init_W(double *R, double *v){

    double sigma_x, sigma_v, mean_x, mean_v;
    int i;

    /* Gaussian number generator  for (R,P) */

    randnums(4*N_bath, ranVector);
    for (i = 0; i< N_bath; i++){
        mean_x = meann[i];
        sigma_x = sig[i];
        // sigma_x = 1.0;
        mean_v = 0;
        sigma_v = sig[i+N_bath];
        // sigma_v = 1.0;
        R[i] = gauss1(mean_x, sigma_x,i);
        v[i] = gauss1(mean_v, sigma_v,i + 2*N_bath);
    }

}


/* For the definition of gamma, dgamma etc, please
   see our JCP paper  */


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

    return x/2.0;


}


void dgamma(double *R){

    int i;

    for (i = 0; i < N_bath; i++)
        dgam[i] = -c[i];

}

void Fb(double *R){

    /* Pure Bath Force Field */

    int i;
    double x;

    for (i= 0; i < N_bath; i++)
        f[i] = mww[i]*R[i];

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

    for (i = 0; i< N_bath; i++)
        f[i] = mww[i]*R[i] + h*c[i];

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
    for (i = 0; i < N_bath; i++)
        dhat[i] /= abs_d;
    /* Note may need to change the sign of dhat here

    /* dhat now is not normalized, and equal d_{10}(R)
       -------------------------------------------------------------- */


}


void (*force[4])(double *);


void integ_step(double *r, double *v, double dt, int Sa){

    /* ********* Velocity Verlet ************************** */

    double y;
    int i;

    y = 0.5*dt*dt;
    for (i = 0; i < N_bath; i++)
        r[i] += dt*v[i] + y*f[i];
    y = 0.5*dt;
    for (i = 0; i < N_bath; i++)
        v[i] +=   y*f[i];
    force[Sa](r);
    for (i = 0; i < N_bath; i++)
        v[i] +=   y*f[i];


}


void bath_para(double eta, double w_max){

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


double U( double *r,double *v, int Sa, double t){

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

int  density(double *x,double *p){


    int l,i,j,SS0,SS1,SS2,SS3, abc,NNjmp =0,nojumpflag = -1,signPdotdhat,adiabat_flag = -1,zcut=0;
    double phase0 = 0.0,phase1 = 0.0,xx,phase3=0,yy;
    double p0,p1,p2,p3,ap0,ap1,ap2,ap3,wtemp = 1.0;
    double dn1,dn2,dn3,bb,sqrbb,sqrbb2,sqrbbb,BETA,GAMMA,pbb1,pbb2,pbb3,ABSZ,ARGZ;
    complex<double> z = 1.0,oldz;

    Dt = TSLICE;
    /* Note - sign for propagating observables  */

    // Initialization of sample


    gauss_init_W(x, p);
    yy = 4.0*(gsl_rng_uniform (rr));
    if (yy < 1.0)
        SS3 = (SS0 = 0);
    else if (yy < 2.0){
        SS0 = 1;
        SS3 = 2;
    }
    else if (yy < 3.0){
        SS0 = 2;
        SS3 = 1;
    }
    else
        SS3 = (SS0 = 3);
    SS[0] = SS0;
    z = 4.0;
    for (l = 0; l < N_bath; l++){
        RR[0][l] = x[l];
        PP[0][l] = p[l];
    }
    SS1 = SS0;

    // ____________________________________________________________________

    adiabat_flag = -1;
    for (l =0; l < N_slice; l++){
        SS0   = SS1;
        // phase0 = 0.0;
        phase0 =  U(RR[0],PP[0],SS0,Dt/2);
        z *= exp(I*phase0);

        //  compute Q1   choos
        //  choose ad or Ua or Ub
        // if Ua or Ub choose kick
        // choose s2 for ad, Ua or Ub
        // note that if ad, s2 = s1
        //  propagate with expLd_{s2}
        // update all prob weights
        // update cummulat for current slice
        // go back to beginimg

        // ****************************************
        //    compute Q1   and choose s1

        dd(RR[0]);
        de = dE(RR[0]);  // note the argument here was x before 5-8-2004
        abc = 0;
        alpha = 0.0;
        bb = de*Dt*abs_d;

        for (Pdotdhat = 0, i = 0; i < N_bath; i++){
            Pdotdhat +=  PP[0][i]*dhat[i];
        }
        alpha = Pdotdhat*abs_d*Dt;
        if (NNjmp >= Ncut)
            adiabat_flag = 1;
        if ( adiabat_flag == 1){
            nojumpflag = 1;
        }

        // 21-8-05 remove 1/p.d exit condition
        if ( 2 > 1){
            signPdotdhat = (Pdotdhat < 0 ? -1 : 1);
            Pdotdhat = fabs(Pdotdhat);
            for (i = 0; i < N_bath; i++)
                Pperp[i] =  PP[0][i] - signPdotdhat*Pdotdhat*dhat[i];
            alpha *= 2.0;
            sina = sin(alpha); cosa = cos(alpha);

            ap0 = fabs(p0 = ( (www[1][SS0][0]() < -7775.0) ? 0.0 :  www[0][SS0][0]()));
            ap1 = fabs(p1 = ( (www[1][SS0][1]() < -7775.0) ? 0.0 :  www[0][SS0][1]()));
            ap2 = fabs(p2 = ( (www[1][SS0][2]() < -7775.0) ? 0.0 :  www[0][SS0][2]()));
            ap3 = fabs(p3 = ( (www[1][SS0][3]() < -7775.0) ? 0.0 :  www[0][SS0][3]()));
            dn2 =  ap0 + ap1 + ap2  + ap3;

            xx = dn2*(gsl_rng_uniform (rr));   // choosing matrix elements
            // printf(" oldz  real %lf imag %lf\t", real(z), imag(z));
            oldz = z;
            SS2 = SS1;
            if (xx < ap0){
                SS1 = 0;
                z *= p0*dn2/ap0;
            }
            else if (xx < ap0 + ap1){
                SS1 = 1;
                z *= p1*dn2/ap1;
            }
            else if (xx < ap0 + ap1 + ap2){
                SS1 = 2;
                z *= p2*dn2/ap2;
            }
            else {
                SS1 = 3;
                z *= p3*dn2/ap3;
            }
            if (SS0 != SS1)
                NNjmp++;
            if (NNjmp > Ncut)
                return 0;

            //  if  ((abs(z) > ppower) || (abs(z) >  4.5*pow(1.5,Dt*l))){
            // 14-09-2006
            if  ((abs(z) > ppower)){
                if (SS0 != SS1)
                    NNjmp--;
                SS1 = SS2;
                z = oldz*www[0][SS0][SS1]();
                // z = oldz;
                // goto jmp;
            }
            if  (www[1][SS0][SS1]() != 9999.0)
                for (i = 0; i < N_bath; i++)
                    PP[0][i] = Pperp[i] + signPdotdhat*www[1][SS0][SS1]()*dhat[i];

        }

        //  printf("zzz %lf %d  %lf\n", l*Dt, NNjmp, abs(z));

        jmp:
        phase0 = U(RR[0],PP[0],SS1,Dt/2);
        z *= exp(I*phase0);
        phi = obs[SS1];


        abszsum0[l]  = real(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        argzsum0[l]  = imag(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        ABSZ = abs(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        ARGZ = arg(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        printf("zzzz %d %d %d %lf %lf\t %lf %lf\t %lf %lf\t %lf %lf\n", l, NNjmp, SS1, real(z), imag(z),  abs(z), arg(z), ABSZ, ARGZ, alpha, de);
        realsum[l][NNjmp] +=  abszsum0[l];
        imagsum[l][NNjmp] +=  argzsum0[l];
        abszsum1[l] += abszsum0[l];
        argzsum1[l] += argzsum0[l];

        hist[l][NNjmp]++;

        phi = obs1[SS1];
        habszsum0[l]  = real(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        hargzsum0[l]  = imag(z*phi(RR[0],PP[0],0)*dens_init[SS3](x,p,1));
        hrealsum[l][NNjmp] +=  habszsum0[l];
        himagsum[l][NNjmp] +=  hargzsum0[l];
        habszsum1[l] += habszsum0[l];
        hargzsum1[l] += hargzsum0[l];

    }



    return 0;

}



int monte(int NN_sample,double *x, double *p){

    int i, j,k,flag,sl;
    double y,l;

    /* This section calls slice algorithm, and prints results  */
    for (i = 0; i < N_slice; i++){
        abszsum1[i] = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;
    }
#pragma omp parallel for num_threads(8)
    {
        for (i = 0; i < NN_sample; i++){
            density(x,p);
            if (((i+1) % Nblock) == 0){
                l  = 1.0/(i+1);
                stream = fopen(datafilename,"a");
                for (k = 0; k < N_slice; k++)
                    if ( ((k+1) % t_strobe) == 0) {
                        for (j = 0; j <=(Ncut+1); j++){
                            fprintf(stream,"%d %lf %d %lf %lf %lf %lf  %lf %lf %lf %lf%.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l, realsum[k][j]*l, imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l );
                            printf("%d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l,realsum[k][j]*l ,imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l);              }
                    }
                fclose(stream);
            }
        }
    }
    return 0;

}

//    Transition Matrices ___________________________________________________________________________-


/* Q1 */



double wwa0_00(){

    double x;

    x = 1 + cosa;
    return x/2.0;

}

double wwa0_01(){

    double x;

    x = -sina;
    return x/2.0;

}

double wwa0_02(){

    double x;

    x = -sina;
    return x/2.0;

}



double wwa0_03(){

    double x;

    x = 1.0 -cosa;
    return x/2.0;

}

double wwa0_10(){

    double x;

    x =  sina;
    return x/2.0;

}

double wwa0_11(){

    double x;

    x = 1.0 + cosa;
    return x/2.0;

}

double wwa0_12(){

    double x;

    x = -1.0 + cosa;
    return x/2.0;

}

double wwa0_13(){

    double x;

    x = -sina;
    return  x/2.0;

}

double wwa0_20(){

    double x;

    x =  sina;
    return x/2.0;

}


double wwa0_21(){

    double x;

    x = -1.0 + cosa ;
    return x/2.0;

}


double wwa0_22(){

    double x;

    x = 1.0 + cosa;
    return x/2.0;

}


double wwa0_23(){

    double x;

    x = -sina ;
    return x/2.0;

}

double wwa0_30(){

    double x;

    x = 1- cosa;
    return x/2.0;

}

double wwa0_31(){

    double x;

    x =  sina;
    return x/2.0;

}

double wwa0_32(){

    double x;

    x =  sina;
    return x/2.0;

}


double wwa0_33(){

    double x;

    x = 1 + cosa;
    return x/2.0;

}

//         W_{a 1}

/* _____________________________________________  */

double wwa1_00(){



    return 9999.0;

}


double wwa1_01(){

    double x;

    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);

}

double wwa1_02(){
    double x;

    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);



}


double wwa1_03(){

    double x;

    x = Pdotdhat*Pdotdhat - 2.0*de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);


}


double wwa1_10(){

    double x;

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

    double x;

    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);


}



double wwa1_20(){

    double x;

    x = Pdotdhat*Pdotdhat + de;

    return sqrt(x);


}

double wwa1_21(){

    double x;


    return 9999.0;

}

double wwa1_22(){

    double x;

    return 9999.0;

}

double wwa1_23(){

    double x;

    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);



}

double wwa1_30(){

    double x;

    x = Pdotdhat*Pdotdhat + 2.0*de;
    return sqrt(x);


}


double wwa1_31(){

    double x;

    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);


}


double wwa1_32(){

    double x;

    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);




}


double wwa1_33(){

    return 9999.0;



}

//   W_{a2}

/* _____________________________________________  */

double wwa2_00(){

    double x;

    x = sina -  sinb2;
    return x;

}


double wwa2_01(){

    double x;

    x = cosb2;
    return x;

}

double wwa2_02(){

    double x;

    x = cosb2;
    return x;

}


double wwa2_03(){

    double x;

    x = sina  + sinb2;
    return  x;

}


double wwa2_10(){

    double x;

    x = cosa;
    return x;

}

double wwa2_11(){

    double x;


    return 0.0;

}

double wwa2_12(){

    double x;

    return 0.0;

}

double wwa2_13(){

    double x;

    x = cosa;
    return x;

}



double wwa2_20(){

    double x;

    x = cosa;
    return x;

}

double wwa2_21(){

    double x;


    return 0.0;

}

double wwa2_22(){

    double x;

    return 0.0;

}

double wwa2_23(){

    double x;

    x = cosa;
    return x;

}

double wwa2_30(){

    double x;

    x = sina +  sinb2;
    return  -x;

}


double wwa2_31(){

    double x;

    x = cosb2;
    return x;

}


double wwa2_32(){

    double x;

    x = cosb2;
    return x;

}


double wwa2_33(){

    double x;

    x = -sina  + sinb2;
    return  x;

}

//   Wb0

double wwb0_00(){

    double x;

    x = cosa;
    return x*x;

}

double wwb0_01(){

    double x;

    x = cosa*sina;
    return x;

}

double wwb0_02(){

    double x;

    x =  cosa*sina;
    return x;

}

double wwb0_03(){

    double x;

    x = sina;
    return 1.0 + x*x;

}

double wwb0_10(){

    double x;

    x = -cosa*sina;
    return x;

}

double wwb0_11(){

    double x;

    x = cosa;
    return x*x;

}

double wwb0_12(){

    double x;

    x = cosa;
    return x*x;

}


double wwb0_13(){

    double x;

    x = cosa*sina;
    return -x;

}


double wwb0_20(){

    double x;

    x = cosa*sina;
    return -x;

}



double wwb0_21(){

    double x;


    x = cosa;
    return x*x;

}

double wwb0_22(){

    double x;

    x = cosa;
    return x*x;

}


double wwb0_23(){

    double x;

    x = cosa*sina;
    return  x;

}


double wwb0_30(){

    double x;

    x = sina;
    return 1.0 + x*x;

}

double wwb0_31(){

    double x;

    x = cosa*sina;
    return -x;

}

double wwb0_32(){

    double x;

    x = -cosa*sina;
    return x;

}


double wwb0_33(){

    double x;

    x = cosa;
    return  x*x;

}



//  W_{b1}


double wwb1_00(){

    double x;

    x = 1.0 - sina*sing1;
    return x;

}


double wwb1_01(){

    double x;

    x = cosg1*sina;
    return x;

}

double wwb1_02(){

    double x;

    x = cosg1*sina;
    return x;

}


double wwb1_03(){

    double x;

    x = 1.0 + sina*sing1;
    return x;

}


double wwb1_10(){

    double x;

    x = cosa*sing1;
    return -x;

}

double wwb1_11(){

    double x;

    x = cosa*cosg1;


    return x;

}

double wwb1_12(){

    double x;

    x = cosa*cosg1;


    return x;


}

double wwb1_13(){

    double x;

    x = cosa*sing1;
    return x;



}



double wwb1_20(){

    double x;

    x = cosa*sing1;
    return -x;

}

double wwb1_21(){

    double x;

    x = cosa*cosg1;


    return x;

}

double wwb1_22(){

    double x;

    x = cosa*cosg1;


    return x;


}

double wwb1_23(){

    double x;

    x = cosa*sing1;
    return x;



}

double wwb1_30(){

    double x;

    x = 1.0 + sina*sing1;
    return x;

}


double wwb1_31(){

    double x;

    x = -cosg1*sina;
    return x;

}

double wwb1_32(){

    double x;

    x = -cosg1*sina;
    return x;

}


double wwb1_33(){

    double x;

    x = 1.0 - sina*sing1;
    return x;

}


//  W_{b2}


double wwb2_00(){

    double x;

    x = 1.0 - sina*sing2;
    return x;

}


double wwb2_01(){

    double x;

    x = cosg2*sina;
    return x;

}

double wwb2_02(){

    double x;

    x = cosg2*sina;
    return x;

}


double wwb2_03(){

    double x;

    x = 1.0 + sina*sing2;
    return x;

}


double wwb2_10(){

    double x;

    x = cosa*sing2;
    return -x;

}

double wwb2_11(){

    double x;

    x = cosa*cosg2;


    return x;

}

double wwb2_12(){

    double x;

    x = cosa*cosg2;


    return x;


}

double wwb2_13(){

    double x;

    x = cosa*sing2;
    return x;



}



double wwb2_20(){

    double x;

    x = cosa*sing2;
    return -x;

}

double wwb2_21(){

    double x;

    x = cosa*cosg2;


    return x;

}

double wwb2_22(){

    double x;

    x = cosa*cosg2;


    return x;


}

double wwb2_23(){

    double x;

    x = cosa*sing2;
    return x;



}

double wwb2_30(){

    double x;

    x = 1.0 + sina*sing2;
    return x;

}


double wwb2_31(){

    double x;

    x = cosg2*sina;
    return -x;

}

double wwb2_32(){

    double x;

    x = cosg2*sina;
    return -x;

}


double wwb2_33(){

    double x;

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




/* ______________________________________   */

/* ************** Observables and initial density Matrices *********  */


/* Definition of initial density matrix element */


double wigner_harm_osc(double *x, double *p){

    double prod,y1,y2;

    int i;

    /* set wigner to one has it drops out of calculation for our initial
    density */

    return 1.0;


}



double dens_init_0(double *x,double *p, int k){

    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);

    return (z*wigner_harm_osc(x,p));


}

double dens_init_1(double *x,double *p, int k){


    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);

    return (z*wigner_harm_osc(x,p));

}

double dens_init_2(double *x,double *p, int k){

    double z;
    double g, gg;


    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);


    return (z*wigner_harm_osc(x,p));

}

double dens_init_3(double *x,double *p, int k){

    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);

    return (z*wigner_harm_osc(x,p));


}

double obs_0(double *x,double *p, int k){


    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 2.0*g/(1 + gg);
    return z;

}

double obs_1(double *x,double *p, int k){


    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = (gg-1)/(1 + gg);
    return z;

}

double obs_2(double *x,double *p, int k){


    double z;
    double g, gg;

    g = G(x); gg = g*g;
    z =  (gg-1)/(1 + gg);
    return z;

}

double obs_3(double *x ,double *p, int k){


    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = -2.0*g/(1 + gg);
    return z;

}

/* These are the matrix elements of the Hamiltonian */

double H_0(double *x,double *p, int k){


    double z;
    double g,gg;

    z = Hb(x,p) - dE(x)/2.0;
    return z;


}



double H_1(double *x,double *p, int k){

    double z;
    double g,gg;
    return 0.0;

}

double H_2(double *x,double *p, int k){

    double z;
    double g, gg;
    return 0.0;


}

double H_3(double *x,double *p, int k){


    double z;
    double g,gg;

    z = Hb(x,p) + dE(x)/2.0;
    return z;

}



int main(int argc, char *argv[]){

    double z,sum1,sum2,z2,sum3,sum5,xx,Dt;
    double tw[4][4];
    double  dt, t1, *R1,  *v,x;
    double w_max, eta,reS, imS,prod,Ndev,yy,T;
    int  i,Sa,j, Nsteps,  N_space,ii,k,init_seed,ppp, N0,stepnum, error_flag= 0;
    int i1, i2;
    t_strobe = 1;   /* Print out the results for every t_strobe slice  */

    gsl_rng_env_setup();

    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);

    // rngstreamnum = 0;
    // nrngstreams = 1;
    // rng_seed = make_sprng_seed();
    // rngstream = init_sprng(rngstreamnum,nrngstreams,rng_seed,SPRNG_DEFAULT);
    /* initialize stream   */
    printf(" Print information about new stream:\n");
    // print_sprng(rngstream);
    fputs("Input datafilename, N_bath, N_slice, Ncut\n timestep, T, init_seed, Nsample\n w_max, eta, beta, delta,power\n ",stderr);
    scanf("%s%d%d%d%lf%lf%d%d%lf%lf%lf%lf%lf", datafilename, &N_bath, &N_slice, &Ncut, &timestep, &T, &init_seed, &Nsample, &w_max, &eta, &beta, &delta, &ppower);
    reS = 0.0; imS = 0.0;
    MAX_tr = 2;
    /* Allocate memory ----------------------------------------------------- */
    R1 = (double *)malloc((N_bath)*sizeof(double));
    v =  (double *)malloc((N_bath)*sizeof(double));
    f =  (double *)malloc((N_bath)*sizeof(double));
    m =  (double *)malloc((N_bath)*sizeof(double));
    w =  (double *)malloc((N_bath)*sizeof(double));
    c =  (double *)malloc((N_bath)*sizeof(double));
    dhat =  (double *)malloc((N_bath)*sizeof(double));
    dgam =  (double *)malloc((N_bath)*sizeof(double));
    b  =  (double *)malloc((N_bath)*sizeof(double));
    Pperp  =  (double *)malloc((N_bath)*sizeof(double));
    SS = (int *)malloc((2)*sizeof(double));
    RR = dmatrix(0,MAX_tr-1,0,N_bath-1);
    PP = dmatrix(0,MAX_tr-1,0,N_bath-1);
    realsum = dmatrix(0,N_slice,0,N_slice+1);
    imagsum = dmatrix(0,N_slice,0,N_slice+1);
    hrealsum = dmatrix(0,N_slice,0,N_slice+1);
    himagsum = dmatrix(0,N_slice,0,N_slice+1);

    hist = imatrix(N_slice,N_slice+1);

    mu =  (double *)malloc((N_bath)*sizeof(double));
    mww = (double *)malloc((N_bath)*sizeof(double));
    dtdtm = (double *)malloc((N_bath)*sizeof(double));
    sig =  (double *)malloc((2*N_bath)*sizeof(double));
    meann = (double *)malloc((2*N_bath)*sizeof(double));
    abszsum0  = (double  *)malloc((N_slice)*sizeof(double));
    abszsum1  = (double  *)malloc((N_slice)*sizeof(double));
    argzsum0  = (double  *)malloc((N_slice)*sizeof(double));
    argzsum1  = (double  *)malloc((N_slice)*sizeof(double));

    habszsum0  = (double  *)malloc((N_slice)*sizeof(double));
    habszsum1  = (double  *)malloc((N_slice)*sizeof(double));
    hargzsum0  = (double  *)malloc((N_slice)*sizeof(double));
    hargzsum1  = (double  *)malloc((N_slice)*sizeof(double));

    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;
    // obs[0] = H_0; obs[1] = H_1; obs[2] = H_2; obs[3] = H_3;
    ddd4 = delta*delta/4.0;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;
    Dt = TSLICE;
    Ndev = 4.0;
    bath_para(eta,w_max);       /* compute system parameters etc */
    /*
    bath corresponds to eq. 54
    for (i = 0; i < N_bath; i++)
       mu[i] = beta*w[i]/2.0;
    for (i = 0; i < N_bath; i++){
       sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
       meann[i] = c[i]/(w[i]*w[i]);
       mww[i] = -m[i]*w[i]*w[i];
       dtdtm[i] = -0.5*timestep*timestep/m[i];
    }

    for (i = 0; i < N_bath; i++){
       sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));
       meann[i+N_bath] = 0.0;
    }
    */
    //  bath corresponds to eq. 53
    for (i = 0; i < N_bath; i++)
        mu[i] = beta*w[i]/2.0;
    for (i = 0; i < N_bath; i++){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
        dtdtm[i] = -0.5*timestep*timestep/m[i];
    }

    for (i = 0; i < N_bath; i++)
        sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));
    force[0] = F1;         /* assign pointers to force fields */
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;
    setwww();
    /*
    alpha = 0.0; BETA = 1.5*alpha; GAMMA = 2.0*alpha;


      for (i1 = 0; i1 < 4 ;i1++)
	for(i2 = 0; i2 < 4; i2++){
           tw[i1][i2] = 0.0;
	   for (j=0;j <3;j++)
	      tw[i1][i2] +=  www[j][i1][i2]();
        }
     for (i1 = 0; i1 < 4 ;i1++){
        printf("\n");
        for(i2 = 0; i2 < 4; i2++)
	    printf("%lf\t", tw[i1][i2]);
        printf("\n");
    }
     printf("OK\n");

    Checked for alpha = beta = gamma = 0 and gives unity as it should

    */


    stream = fopen(datafilename,"w");
    fprintf(stream,"%s\n w_max %lf eta %lf beta %lf delta %lf killz %lf N_bath %d N_slice %d\n", argv[0], w_max, eta, beta, delta, ppower, N_bath, N_slice);
    fprintf(stream,"Nens\t time\t j\t O_j\t O_tot\t En_j\t En_tot\n");
    fclose(stream);
    for (i = 0; i < N_slice; i++)
        for (j = 0; j <= N_slice; j++){
            realsum[i][j] = 0.0;
            imagsum[i][j] = 0.0;
            hrealsum[i][j] = 0.0;
            himagsum[i][j] = 0.0;


            hist[i][j] = 0;
        }
    monte(Nsample,R1,v);
    stream = fopen(datafilename,"a");
    fprintf(stream,"dt %lf T %lf Nsample\n", timestep, T, Nsample);
    for (i = 0; i < N_slice; i++)
        if (((i+1)% t_strobe) == 0)
            for (j =0; j <= (Ncut+1);j++){
                printf("%lf   %lf  %lf  %lf  %lf  %lf  %lf    %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample,  (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample,  1.0*hist[i][j]/Nsample);
                fprintf(stream,"%lf   %lf  %lf  %lf  %lf  %lf  %lf %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample, (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample, 1.0*hist[i][j]/Nsample);


            }
    fclose(stream);

    free(R1);  free(v); free(f); free(b);
    free(m), free(c); free(w); free(dhat); free(dgam);
    free(SS); free(jump); free(Pperp);
    free_dmatrix(RR,0,MAX_tr-1,0,N_bath-1);
    free_dmatrix(PP,0,MAX_tr-1,0,N_bath-1);
    free_dmatrix(realsum,0,N_slice,0,N_slice+1);
    free_dmatrix(imagsum,0,N_slice,0,N_slice+1);
    free_dmatrix(hrealsum,0,N_slice,0,N_slice+1);
    free_dmatrix(himagsum,0,N_slice,0,N_slice+1);
    free_imatrix(hist,N_slice,N_slice+1);
    free(abszsum0);  free(argzsum0);    free(abszsum1);  free(argzsum1);
    free(habszsum0);  free(hargzsum0);    free(habszsum1);free(hargzsum1);
    free(mu);free(sig); free(mww);free(dtdtm); free(meann);

    return 0;

}