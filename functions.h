#ifndef FUNCTIONS
#define FUNCTIONS

double gam(double *R);
double Hb(double *R, double *P);
void Fb(double *R);
void F1(double *R);
void F2(double *R);
double dE(double *R);
double G(double *R);
void dd(double*R);
void integ_step(double *r, double *v, double dt, int Sa);
void bath_para(double eta, double w_max);
double U( double *r,double *v, int Sa, double t);


#endif
