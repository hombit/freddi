#ifndef _SPECTRUM_HPP
#define _SPECTRUM_HPP



#include <cmath>
#include <vector>

#include "gsl_const_cgsm.h"


double Luminosity(const std::vector<double> &R, const std::vector<double> &T, double min_nu, double max_nu, int Nnu);

double I_lambda( const std::vector<double> &R, const std::vector<double> &T, double lambda );

double T_GR( const double r1, const double ak, const double Mx, const double Mdot, const double r_in  );


#endif // _SPECTRUM_HPP
