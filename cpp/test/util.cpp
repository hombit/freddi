#include <cmath>
#include <vector>

#include <util.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_util

#include <boost/test/unit_test.hpp>


std::vector<double> get_x(size_t N) {
	const double x_min = 2 * M_PI;
	const double x_max = 3 * M_PI;

	std::vector<double> x(N);
	for (size_t i = 0; i < N; i++) {
		x[i] = x_min * std::pow(x_max / x_min, i / (N - 1.));
	}
	return x;
}

std::vector<double> get_y(const std::vector<double>& x) {
	std::vector<double> y(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		y[i] = std::sin(x[i]);
	}
	return y;
}

BOOST_AUTO_TEST_CASE(test_trapz) {
	const size_t N = 11;
	const auto x = get_x(N);
	const auto y = get_y(x);

	const auto result = trapz(x, y, 0, N-1);
	// Analytical
	BOOST_CHECK_CLOSE_FRACTION(result, 2, 1e-2);
	// scipy.integrate.trapz
	BOOST_CHECK_CLOSE_FRACTION(result, 1.9829331321624648, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_simps_odd) {
	const size_t N = 11;
	const auto x = get_x(N);
	const auto y = get_y(x);

	const auto result = simps(x, y, 0, N-1);
	// Analytical
	BOOST_CHECK_CLOSE_FRACTION(result, 2, 1e-4);
	// scipy.integrate.simps
	BOOST_CHECK_CLOSE_FRACTION(result, 2.0001763497127496, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_simps_2) {
	const std::vector<double> x({0, 1});
	const std::vector<double> y({2, 4});
	const auto result = simps(x, y, 0, 1);
	BOOST_CHECK_EQUAL(result, 3);
}

BOOST_AUTO_TEST_CASE(test_simps_even) {
	const size_t N = 12;
	const auto x = get_x(N);
	const auto y = get_y(x);

	const auto result = simps(x, y, 0, N-1);
	// Analytical
	BOOST_CHECK_CLOSE_FRACTION(result, 2, 1e-5);
	// scipy.integrate.simps(even='right')
	BOOST_CHECK_CLOSE_FRACTION(result, 1.9999976227091623, 1e-12);
}
