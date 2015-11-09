#ifndef CONVOLUTION_H
#define CONVOLUTION_H

namespace jun
{

void convolution(double x[],
                 double y[],
                 double cx[],
                 double cy[],
                 double xmin,
                 double xmax,
                 int n,
                 int n_out,
                 double sigma);

}

#endif
