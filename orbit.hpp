#ifndef _ORBIT_HPP
#define _ORBIT_HPP


#include <cmath>
#include <gsl/gsl_const_cgsm.h>


double r_out_func(const double Mx, const double Mopt, const double P);

double r_ISCO(const double kerr);

double r_in_func(double Mx, double kerr);

double efficiency_of_accretion(const double kerr);


#endif // _ORBIT_HPP
