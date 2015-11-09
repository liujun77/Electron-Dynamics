#include <complex>
#include "libAmosBessel/libAmosBessel.h"
#include <iostream>

using namespace std;

namespace jun{

complex<double> bessel_J_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('J',z,n,1,0,&f);
    return f;
}

complex<double> bessel_Y_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('Y',z,n,1,0,&f);
    return f;
}

complex<double> bessel_I_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('I',z,n,1,0,&f);
    return f;
}

complex<double> bessel_K_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('K',z,n,1,0,&f);
    return f;
}

complex<double> bessel_j_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('j',z,n,1,0,&f);
    return f;
}

complex<double> bessel_y_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('y',z,n,1,0,&f);
    return f;
}

complex<double> bessel_h1_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('o',z,n,1,0,&f);
    return f;
}

complex<double> bessel_h2_n(double n, complex<double> z){
    int i;
    complex<double> f;
    i=AmosBessel('t',z,n,1,0,&f);
    return f;
}

}
