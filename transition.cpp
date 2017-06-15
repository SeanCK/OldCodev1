#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "variables.cpp"


using namespace std;

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

    /*   note dE = E1 - E0 > 0 , w_10 = dE , corresponding to the 1 0 pair
     of surfaces, equivalent to Sa = 2    */


double alpha,BETA1,BETA2, GAMMA1, GAMMA2,*b;



double sina, cosa, cosb1,cosb2, sinb1,sinb2,cosg1,cosg2,sing1,sing2, ppower,de;


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
