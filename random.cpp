#include   <stdio.h>
#include   <stdio.h>
#include   <stdlib.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "variables.cpp"

using namespace std;

#include <gsl/gsl_rng.h>

#define PI 3.141592653589793


/* Random number stuff ***************************************************************  */


double  gauss1(double mean_x, double sigma_x,int i){

    double x1,y1,y2,y3;

    /* The tree here is to avoid a minor instability ***************************       */

    y1 = ranVector[i];

    if (fabs(y1) < 1.0e-200){ // Why are there 4 iterations? Do you want a number < 1.0e-200? If >, leaves the loop
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
    return (mean_x = sigma_x*x1);
}


void randnums(int rand_dim, double *rand_vec){ //No returned value

    int i;
    for (i = 0; i < rand_dim; i++){
        rand_vec[i] = gsl_rng_uniform (rr);
    }
}


void gauss_init_W(double *R, double *v){ //No returned value

    double sigma_x, sigma_v, mean_x, mean_v;
    int i;

    /* Gaussian number generator  for (R,P) */

    randnums(4*N_bath, ranVector);
    for (i = 0; i < N_bath; i++){
        mean_x = meann[i];
        sigma_x = sig[i];
        // sigma_x = 1.0;
        mean_v = 0;
        sigma_v = sig[i+N_bath];
        // sigma_v = 1.0;
        R[i] = gauss1(mean_x, sigma_x, i);
        v[i] = gauss1(mean_v, sigma_v, i + 2*N_bath);
    }
}

int  density(double *x,double *p){


    int l,i,j,SS0,SS1,SS2,SS3, abc, NNjmp =0, nojumpflag = -1, signPdotdhat, adiabat_flag = -1, zcut=0;
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