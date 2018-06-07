#ifndef _CONSTANTS_HPP
#define _CONSTANTS_HPP


const double DAY = 86400.;
const double Angstrem = 1e-8;
const double keV = 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
const double Jy = 1e-23;
const double solar_radius = 6.955e10;
const double kpc = 1000. * GSL_CONST_CGSM_PARSEC;

namespace photometry_bands{
// Allen's Astrophysical Quantities (4th ed.)
const double lambdaU = 3600. * Angstrem;
const double irr0U = 4.22e-9 / Angstrem;
const double lambdaB = 4400. * Angstrem;
const double irr0B = 6.4e-9 / Angstrem;
const double lambdaV = 5500. * Angstrem;
const double irr0V = 3.750e-9 / Angstrem;
const double lambdaR = 7100 * Angstrem;
const double irr0R = 1.75e-9 / Angstrem;
const double lambdaI = 9700 * Angstrem;
const double irr0I = 0.84e-9 / Angstrem;
// Campins et al., 1985, AJ, 90, 896
const double lambdaJ = 12600 * Angstrem;
const double irr0J = 1600 * Jy *  GSL_CONST_CGSM_SPEED_OF_LIGHT / (lambdaJ*lambdaJ); 
} // photometry_bands

#endif // _CONSTANTS_HPP
