// *****************************************************************************
// *** Constants for Physics                                                 ***
// *** Jun Liu 2014-02-09                                                    ***
// *****************************************************************************

#include <complex>

namespace jun{

extern const double PI=3.141592653589793238462643;
extern const double AU_EV=27.2113834;        //1 a.u. of energy in eV
extern const double A0_NM=0.0529177208;      //1 Bohr radius in nm
extern const double C_AU=137.03599971;       //1/alfa=speed of light in a.u.
extern const double Boltzman_k=8.6173324e-5/27.211; //Boltzman const in au/K
extern const double EC=1.602176565e-19;      //Elementary charge in C

//Factorial generation
long *FACTORIAL;

int CONST_factorial_generator(unsigned n){
    unsigned i;
    FACTORIAL=new long[n];
    FACTORIAL[0]=1L;
    for(i=1;i<=n;++i){
        FACTORIAL[i]=FACTORIAL[i-1]*i;
    }
    return 1;
}

double *DFACTORIAL;

int CONST_dfactorial_generator(unsigned n){
    unsigned i;
    DFACTORIAL=new double[n];
    DFACTORIAL[0]=1.;
    for(i=1;i<=n;++i){
        DFACTORIAL[i]=DFACTORIAL[i-1]*i;
    }
    return 1;
}

std::complex<double> I(0,1);

}
