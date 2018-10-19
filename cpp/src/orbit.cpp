#include <cmath>

#include "gsl_const_cgsm.h"

#include "orbit.hpp"


double rISCORg(double kerr) { // From «Black Hole Accretion Disks», A.44 (p. 530)
	const double Z1 = 1. + std::cbrt((1. - kerr * kerr)) * (std::cbrt((1. + kerr)) + std::cbrt((1. - kerr)));
	double Z2 = std::sqrt(3. * kerr * kerr + Z1 * Z1);
	return 3. + Z2 - std::sqrt((3. - Z1) * (3. + Z1 + 2. * Z2));
}

double efficiency_of_accretion(const double kerr) {
	return 1. - std::sqrt(1. - 2. / 3. / rISCORg(kerr));
}
