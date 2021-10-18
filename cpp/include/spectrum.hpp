#ifndef _SPECTRUM_HPP
#define _SPECTRUM_HPP



#include <cmath>
#include <limits>
#include <vector>

#include "gsl_const_cgsm.h"
#include "util.hpp"


namespace Spectrum {
constexpr const double double_h_over_c2 = 2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
constexpr const double h_over_kB = GSL_CONST_CGSM_PLANCKS_CONSTANT_H / GSL_CONST_CGSM_BOLTZMANN;
constexpr const double double_h_c2 = 2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT;
constexpr const double ch_over_kB = GSL_CONST_CGSM_SPEED_OF_LIGHT * h_over_kB;
double Planck_nu(double T, double nu);
double Planck_lambda(double T, double lambda);

double Planck_nu1_nu2(double T, double nu1, double nu2, double tol=std::sqrt(std::numeric_limits<double>::epsilon()));

double T_GR(double r1, double ak, double Mx, double Mdot);
} // namespace Spectrum


#endif // _SPECTRUM_HPP
