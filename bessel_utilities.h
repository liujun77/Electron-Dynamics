// *****************************************************************************
// *** bessel_utilities.h                                                    ***
// *** advanced bessel functions                                             ***
// *****************************************************************************

#ifndef BESSEL_UTILITIES_H
#define BESSEL_UTILITIES_H

namespace jun{

//besselh1(l,ix)int l, real x
std::complex<double> sph_hankel_i(int l, double x);

//besselh1'(l,ix)int l, real x
std::complex<double> sph_hankel_ip(int l, double x);

//h1'(ix)/h1(ix)
std::complex<double> sph_hankel_ip_h(int l, double x);

double sph_besselp(int l, double x);

//jp/j
double sph_bessel_jp_j(int l, double x);

}

#endif
