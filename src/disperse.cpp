#include <iostream>


#include "disperse.h"

double Dispersal::dist_exponential(xorshift64& rand, double sigma)
{
    double param = 1.0/sigma;
    return rand_exp(rand, param);
}

double Dispersal::dist_triangular(xorshift64& rand, double sigma)
{
    //where xmin = 0 and xmax = c
    double param = 2.0*sigma;
    double u = rand.get_double52();
    return param*sqrt(u);
}

double Dispersal::dist_halfNormal(xorshift64& rand, double sigma)
{
    double param = sigma * sqrt(2.0);
    return rand_abs_normal(rand, 0.0, param);
}

double Dispersal::dist_rayleigh(xorshift64& rand, double sigma)
{
    double param = sigma;
    return param * sqrt(2.0 * rand_exp(rand));
}




