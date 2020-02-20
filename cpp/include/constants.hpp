#ifndef _CONSTANTS_HPP
#define _CONSTANTS_HPP

#include <cmath>

#include "gsl_const_cgsm.h"

constexpr const double FOUR_M_PI = 4.0 * M_PI;

constexpr const double DAY = 86400.;
constexpr const double Angstrem = 1e-8;
constexpr const double keV = 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
constexpr const double Jy = 1e-23;
constexpr const double solar_radius = 6.955e10;
constexpr const double kpc = 1000. * GSL_CONST_CGSM_PARSEC;

// Allen's Astrophysical Quantities (4th ed.)
constexpr const double lambdaU = 3600. * Angstrem;
constexpr const double irr0U = 4.22e-9 / Angstrem;
constexpr const double lambdaB = 4400. * Angstrem;
constexpr const double irr0B = 6.4e-9 / Angstrem;
constexpr const double lambdaV = 5500. * Angstrem;
constexpr const double irr0V = 3.750e-9 / Angstrem;
constexpr const double lambdaR = 7100 * Angstrem;
constexpr const double irr0R = 1.75e-9 / Angstrem;
constexpr const double lambdaI = 9700 * Angstrem;
constexpr const double irr0I = 0.84e-9 / Angstrem;
// Campins et al., 1985, AJ, 90, 896
constexpr const double lambdaJ = 12600 * Angstrem;
constexpr const double irr0J = 1600 * Jy *  GSL_CONST_CGSM_SPEED_OF_LIGHT / (lambdaJ*lambdaJ);

#endif // _CONSTANTS_HPP
