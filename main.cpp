#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "dmatrix.c"
#include   "imatrix.c"
#include   "variables.cpp"
#include   "random.cpp"
#include   "transition.cpp"

using namespace std;
/* #define SIMPLE_SPRNG		 simple interface                         */
// #include "sprng.h"
/* SPRNG header file     see NSCA random number generator   library   *** */


/* For the definition of gamma, dgamma etc, please
   see our JCP paper  */


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

// Monte Carlo Sampling -------------------------------------------------

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

/* ______________________________________   */

/* ************** Observables and initial density Matrices *********  */


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