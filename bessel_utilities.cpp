// *****************************************************************************
// *** bessel_utilities.cpp                                                  ***
// *****************************************************************************

#include <cmath>
#include <complex>
#include "const.h"
#include <boost/math/special_functions/bessel.hpp>

using boost::math::sph_bessel;
using namespace std;
using jun::PI;
using jun::I;
using jun::DFACTORIAL;

namespace jun{

complex<double> sph_hankel_i(int l, double x){ //besselh1(l,ix)int l, real x
    complex<double> factor=pow(I,l);
    int m,s;
    double sum=1,p;
    for(m=1; m<=l; m++){
        p=1.;
        for(s=1; s<=m; s++)
            p=p*2*x;
    sum+=DFACTORIAL[l+m]/DFACTORIAL[m]/DFACTORIAL[l-m]/p;
    }
    l=l%2;
    if(l==0)
        factor=-1.*factor;
    return sum*factor*exp(-x)/x;
}

complex<double> sph_hankel_ip(int l, double x){ //besselh1'(l,ix)int l, real x
    if(l==0)
        return -(1/x/x+1/x)*exp(-x)*I;
    else
        return -0.5*sph_hankel_i(l+1,x)
               +0.5*sph_hankel_i(l-1,x)
               -0.5/I*sph_hankel_i(l,x)/x;
}

complex<double> sph_hankel_ip_h(int l, double x){ //h1'(ix)/h1(ix)
    int m,s;
    long double sum1,sum2,p;
    sum1=0;
    sum2=0;
    for(m=0;m<=l;m++){
        p=DFACTORIAL[l+m]/DFACTORIAL[m]/DFACTORIAL[l-m];
        if(x>0.5)
            for(s=1;s<=m;s++)
                p=p/x/2;
        else
            for(s=1;s<=l-m;s++)
                p=p*2*x;
        sum1+=(m+1)/x*p;
        sum2+=p;
    }
    return I*double(1+sum1/sum2);
}

double sph_besselp(int l, double x){
    if(l==0)
        return -sph_bessel(1,x);
    else
        return sph_bessel(l-1,x)-((l+1)/x)*sph_bessel(l,x);
}

double sph_bessel_jp_j(int l, double x){  //jp/j
    long double jp,j,lx;
    lx=1L*x;
    j=sph_bessel(l,lx);
    if(l==0)
        jp=-sph_bessel(1,x);
    else
        jp=sph_bessel(l-1,x)-((l+1)/x)*j;
  return jp/j;
}

}
