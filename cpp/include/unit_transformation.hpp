#ifndef FREDDI_UNIT_TRANSFOMATION_H
#define FREDDI_UNIT_TRANSFOMATION_H

#include "gsl_const_cgsm.h"
#include "constants.hpp"


constexpr inline double gramToSun(const double mass_gram) { return mass_gram / GSL_CONST_CGSM_SOLAR_MASS; }
constexpr inline double sunToGram(const double mass_sun) { return mass_sun * GSL_CONST_CGSM_SOLAR_MASS; }


constexpr inline double cmToRg(const double length_cm, const double mass_gram) {
	return length_cm * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram);
}
constexpr inline double rgToCm(const double length_rg, const double mass_gram) {
	return length_rg * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
}
constexpr inline double kpcToCm(const double length_kpc) {
	return length_kpc * 1000. * GSL_CONST_CGSM_PARSEC;
}
constexpr inline double cmToKpc(const double length_cm) {
	return length_cm / (1000. * GSL_CONST_CGSM_PARSEC);
}
inline double sunToCm(const double length_solar_radius) {
	return length_solar_radius * solar_radius;
}
constexpr inline double cmToSun(const double length_cm) {
	return length_cm / solar_radius;
}
constexpr inline double angstromToCm(const double length_angstrom) {
	return length_angstrom * 1e-8;
}
constexpr inline double cmToAngstrom(const double length_cm) {
	return length_cm * 1e8;
}


constexpr inline double sToDay(double time_s) { return time_s / 86400.; }
constexpr inline double dayToS(double time_day) { return time_day * 86400.; }


constexpr inline double kevToHertz(const double energy_keV) {
	return energy_keV * 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
}
constexpr inline double hertzToKev(const double nu_hertz) {
	return nu_hertz * GSL_CONST_CGSM_PLANCKS_CONSTANT_H / (1000. * GSL_CONST_CGSM_ELECTRON_VOLT);
}


constexpr inline double kevToK(const double temperature_kev) {
	return temperature_kev * 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_BOLTZMANN;
}
constexpr inline double kToKev(const double temperature_k) {
	return temperature_k * GSL_CONST_CGSM_BOLTZMANN / (1000. * GSL_CONST_CGSM_ELECTRON_VOLT);
}

#endif //FREDDI_UNIT_TRANSFOMATION_H
