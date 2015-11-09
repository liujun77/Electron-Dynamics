// *****************************************************************************
// *** Convolution                                                           ***
// *** Jun Liu 2015-01-26                                                    ***
// *****************************************************************************

#include <cmath>
#include "const.h"

namespace jun
{

double con_gaussian(double x, double sigma, double mu)
{
    return 1./sigma/sqrt(2*jun::PI)*exp(-0.5*(x-mu)*(x-mu)/sigma/sigma);
}

void convolution(double x[],
                 double y[],
                 double cx[],
                 double cy[],
                 double xmin,
                 double xmax,
                 int n,
                 int n_out,
                 double sigma)
{
    double dx=(xmax-xmin)/n_out;
    int i;
    for(i=0;i<n_out;++i){
        cx[i]=xmin+dx*i;
        cy[i]=0;
    }
    for(i=0;i<n_out;++i){
        for(int j=0;j<n;++j)
            cy[i]=cy[i]+y[j]*con_gaussian(cx[i],sigma,x[j]);
    }
}

}
