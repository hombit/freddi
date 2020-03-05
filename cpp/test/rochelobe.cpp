#include <boost/numeric/odeint.hpp>

#include <rochelobe.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testRocheLobe

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace odeint = boost::numeric::odeint;


BOOST_AUTO_TEST_CASE(testRochePolarPotential_derivatives) {
	const std::vector<double> mass_ratio = {1e-6, 0.1, 1.0, 10., 1e6};
	const std::vector<double> lambda_ = {-1.0, 0.5, 0.0, 0.3, 1.0};
	const std::vector<double> nu_ = {-1.0, 0.3, 0.0, 0.5, 1.0};
	const double r1 = 0.1;
	const double r2 = 0.7;
	const double tol = 1e-4;

	auto stepper = odeint::runge_kutta_cash_karp54<double>();
	double integral = 0.;

	for (double q : mass_ratio) {
		for (double lambda : lambda_) {
			for (double nu: nu_) {
				const RochePotential rp(q);

				integral = 0;
				integrate_adaptive(
						odeint::make_controlled(0, tol, stepper),
						[&](const double &y, double &dydx, double x) {
							dydx = rp.domega_dr(x, lambda, nu);
						},
						integral, r1, r2, tol * (r2 - r1)
				);
				BOOST_CHECK_CLOSE(integral, rp.omega(r2, lambda, nu) - rp.omega(r1, lambda, nu), 1e3 * tol);

				integral = 0;
				integrate_adaptive(
						odeint::make_controlled(0, tol, stepper),
						[&](const double &y, double &dydx, double x) {
							dydx = rp.d2omega_dr2(x, lambda, nu);
						},
						integral, r1, r2, tol * (r2 - r1)
				);
				BOOST_CHECK_CLOSE(integral, rp.domega_dr(r2, lambda, nu) - rp.domega_dr(r1, lambda, nu), 1e3 * tol);

				integral = 0;
				integrate_adaptive(
						odeint::make_controlled(0, tol, stepper),
						[&](const double &y, double &dydx, double x) {
							dydx = rp.d3omega_dr3(x, lambda);
						},
						integral, r1, r2, tol * (r2 - r1)
				);
				BOOST_CHECK_CLOSE(integral, rp.d2omega_dr2(r2, lambda, nu) - rp.d2omega_dr2(r1, lambda, nu), 1e3 * tol);
			}
		}
	}
}


BOOST_AUTO_TEST_CASE(testCriticalRocheLobe_lagrangian_point) {
	const std::vector<double> mass_ratio = {1e-3, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3};
	const std::vector<double> expected = {0.06769101046002532, 0.14147491318312078, 0.2824874128854856, 0.5, 0.7175125871145017, 0.8585250868168659, 0.9323089895399626};

	for (size_t i = 0; i < mass_ratio.size(); ++i) {
		CriticalRocheLobe crl(mass_ratio[i]);
		BOOST_CHECK_CLOSE(crl.lagrangian_point, expected[i], 1e-6);
	}
}


BOOST_AUTO_TEST_CASE(testCriticalRocheLobe_polar_radius) {
	const std::vector<double> mass_ratio = {1e-3, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3};
	const std::vector<double> expected = {0.04360976056, 0.09230376785, 0.1898767883, 0.3561308722, 0.5344988591, 0.6286781532, 0.6576154360};

	for (size_t i = 0; i < mass_ratio.size(); ++i) {
		CriticalRocheLobe crl(mass_ratio[i]);
		BOOST_CHECK_CLOSE(crl.polar_radius, expected[i], 1e-5);
	}
}
