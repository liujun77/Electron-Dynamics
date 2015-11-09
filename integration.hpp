// *****************************************************************************
// *** integration.h                                                         ***
// *** Jun Liu                                                               ***
// *****************************************************************************

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <cmath>

namespace jun{

template <typename T1, typename T2>

inline T1 integ(T2 a, T2 b, T1 f(T2 x), int n){
    if(n<7)
        n=7;
    int i;
    T1 sum=0;
    T2 h=(b-a)/n;
    if(a==b)
        return 0;
    sum=(17*(f(a)+f(b))+59*(f(a+h)+f(b-h))
        +43*(f(a+2*h)+f(b-2*h))+49*(f(a+3*h)+f(b-3*h)))/48;
    for(i=4;i<=n-4;++i)
        sum+=f(a+i*h);
    return sum*h;
}

template <typename T1, typename T2, typename T3>

inline T1 integ(T2 a, T2 b, T1 f(T2, T3), T3 ind, int n){
    if(n<7)
        n=7;
    int i;
    T1 sum=0;
    T2 h=(b-a)/n;
    if(a==b)
        return 0;
    sum=(17*(f(a,ind)+f(b,ind))+59*(f(a+h,ind)+f(b-h,ind))
        +43*(f(a+2*h,ind)+f(b-2*h,ind))+49*(f(a+3*h,ind)+f(b-3*h,ind)))/48;
    for(i=4;i<=n-4;++i){
        sum+=f(a+i*h,ind);
#ifdef DEBUG_INTEGRATION_H
        if(std::isnan(sum)){
            std::cout<<"i="<<i
                <<"ind="<<ind<<std::endl;
            std::exit(0);
        }
#endif
    }
    return sum*h;
}

}

#endif
