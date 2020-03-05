#ifndef _ORBIT_HPP
#define _ORBIT_HPP

#include <cmath>

#include <unit_transformation.hpp>
#include <util.hpp>


double r_kerrISCORg(double kerr);

double efficiency_of_accretion(double kerr);


class BlackHoleFunctions {
public:
	static inline double r_kerrISCO(const double Mx, const double kerr) { return rgToCm(r_kerrISCORg(kerr), Mx); }
};


class BinaryFunctions {
public:
	static double rocheLobeVolumeRadiusSemiaxis(double mass_ratio);
	static inline double semiaxis(const double total_mass, const double period) {
		return std::cbrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * total_mass * m::pow<2>(period) / (4. * m::pow<2>(M_PI)));
	}
	static inline double semiaxis(const double mass1, const double mass2, const double period) {
		return semiaxis(mass1 + mass2, period);
	}
	static inline double rocheLobeVolumeRadius(const double mass1, const double mass2, const double period) {
		return rocheLobeVolumeRadiusSemiaxis(mass1 / mass2) * semiaxis(mass1, mass2, period);
	}
};


#endif // _ORBIT_HPP
