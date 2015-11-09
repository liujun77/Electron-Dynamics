#ifndef CONST_H
#define CONST_H

// *****************************************************************************
// *** Constants for Physics                                                 ***
// *** Jun Liu 2014-02-09                                                    ***
// *****************************************************************************

namespace jun{

extern const double PI;//=3.141592653589793238462643,
extern const double AU_EV;//=27.2113834
                    //1 a.u. of energy in eV
extern const double A0_NM;//=0.0529177208
                    //1 Bohr radius in nm
extern const double C_AU;//=137.03599971
                    //1/alfa=speed of light in a.u.
extern const double Boltzman_k;//=8.6173324e-5
                    //Boltzman const in eV/K
extern const double EC;//=1.602176565e-19
                    //Elementary charge in C

//Factorial generation
extern long *FACTORIAL;

int CONST_factorial_generator(unsigned n);

extern double *DFACTORIAL;

int CONST_dfactorial_generator(unsigned n);

#ifdef _GLIBCXX_COMPLEX

extern std::complex<double> I;

#endif

}

#endif
