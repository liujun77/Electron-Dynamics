#ifndef BESSEL_H
#define BESSEL_H

// *****************************************************************************
// *** Bessel functions                                                      ***
// *** Jun Liu 2014-09-16                                                    ***
// *****************************************************************************

#include <complex>

namespace jun{

std::complex<double> bessel_J_n(double n, std::complex<double> z);
std::complex<double> bessel_Y_n(double n, std::complex<double> z);
std::complex<double> bessel_I_n(double n, std::complex<double> z);
std::complex<double> bessel_K_n(double n, std::complex<double> z);
std::complex<double> bessel_j_n(double n, std::complex<double> z);
std::complex<double> bessel_y_n(double n, std::complex<double> z);
std::complex<double> bessel_h1_n(double n, std::complex<double> z);
std::complex<double> bessel_h2_n(double n, std::complex<double> z);

}

#endif
