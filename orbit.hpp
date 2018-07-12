#ifndef _ORBIT_HPP
#define _ORBIT_HPP


#include <cmath>

#include "gsl_const_cgsm.h"


namespace disk_orbit {

double r_out_func(double Mx, double Mopt, double P);

double r_ISCO(double kerr);

double r_in_func(double Mx, double kerr);

double efficiency_of_accretion(double kerr);

} // disk_orbit


#endif // _ORBIT_HPP
