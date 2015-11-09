#ifndef LEGENDRE_RULE_FAST_H
#define LEGENDRE_RULE_FAST_H

void legendre_compute_glr ( int n, double x[], double w[] );

void legendre_compute_glr0 ( int n, double *p, double *pp );

void legendre_compute_glr1 ( int n, double *roots, double *ders );

void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );

void rescale ( double a, double b, int n, double x[], double w[] );

double rk2_leg ( double t, double tn, double x, int n );

double ts_mult ( double *u, double h, int n );

#endif
