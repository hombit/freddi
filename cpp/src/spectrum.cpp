#include "orbit.hpp"
#include "spectrum.hpp"

#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;


namespace Spectrum {
double Planck_nu(const double T, const double nu) {
	if (!(T > 0.)) {  // catches NaN
		return 0.;
	}
	return double_h_over_c2 * m::pow<3>(nu) / (std::exp(h_over_kB * nu / T) - 1.);
}


double Planck_lambda(const double T, const double lambda) {
	if (!(T > 0.)) {  // catches NaN
		return 0.;
	}
	return double_h_c2 / m::pow<5>(lambda) / (std::exp( ch_over_kB / (lambda * T) ) - 1.);
}


double Planck_nu1_nu2(const double T, const double nu1, const double nu2, const double tol){
	auto stepper = odeint::runge_kutta_cash_karp54<double>();
	double integral = 0.;
	integrate_adaptive(
			odeint::make_controlled(0, tol, stepper),
			[T](const double &y, double &dydx, double x) {
				dydx = Planck_nu(T, x);
			},
			integral, nu1, nu2, 1e-2 * (nu2 - nu1)
	);
	return integral;
}


// Code by Galina Lipunova:
/* General Relativity effects are included in the structure of the disk
   (Page & Thorne 1974; Riffert & Herold 1995). metric = "GR"
*/
double T_GR(const double r1, const double ak, const double Mx, const double Mdot){
	const double GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx;
	const double rg = GM  / m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
	const double x = std::sqrt(r1 / rg);
	const double x0 = std::sqrt(r_kerrISCORg(ak));

	if (x < x0) {
		return 0.;
	}
	
	const double x1 = 2. * std::cos ((std::acos(ak)-M_PI)/3.);
	const double x2 = 2. * std::cos ((std::acos(ak)+M_PI)/3.);
	const double x3 = -2. * std::cos (std::acos(ak)/3.);
	const double a = 3. * (x1-ak)*(x1-ak) * std::log((x-x1)/(x0-x1))/x1/(x1-x2)/(x1-x3);
	const double b = 3. * (x2-ak)*(x2-ak) * std::log((x-x2)/(x0-x2))/x2/(x2-x1)/(x2-x3);
	const double c = 3. * (x3-ak)*(x3-ak) * std::log((x-x3)/(x0-x3))/x3/(x3-x1)/(x3-x2);

	return( std::pow(
			(3.*Mdot * m::pow<6>(GSL_CONST_CGSM_SPEED_OF_LIGHT) / (8.*M_PI * m::pow<2>(GM))) *
			(x - x0 - 1.5 * ak * std::log(x/x0) - a - b -c) / ( m::pow<4>(x)*(m::pow<3>(x) - 3.*x + 2. * ak) ) /
			GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT,
			0.25)

	);
}
} // namespace Spectrum
