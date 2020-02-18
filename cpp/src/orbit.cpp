#include <cmath>

#include "gsl_const_cgsm.h"

#include "orbit.hpp"
#include "util.hpp"


double r_kerrISCORg(double kerr) { // From «Black Hole Accretion Disks», A.44 (p. 530)
	const double Z1 = 1. + std::cbrt((1. - m::pow<2>(kerr))) * (std::cbrt((1. + kerr)) + std::cbrt((1. - kerr)));
	double Z2 = std::sqrt(3. * m::pow<2>(kerr) + m::pow<2>(Z1));
	return 3. + Z2 - std::sqrt((3. - Z1) * (3. + Z1 + 2. * Z2));
}

double efficiency_of_accretion(const double kerr) {
	return 1. - std::sqrt(1. - 2. / 3. / r_kerrISCORg(kerr));
}

double BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(const double mass_ratio) {
	// Eggleton, P. P. 1983, ApJ, 268, 368
	const double q = std::cbrt(mass_ratio);
	return 0.49 * m::pow<2>(q) / (0.6 * m::pow<2>(q) + std::log(1. + q));
}
