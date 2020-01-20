#include <boost/math/special_functions/pow.hpp>

#include <constants.hpp>
#include <star.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testStar

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>


namespace m = boost::math;


BOOST_AUTO_TEST_CASE(testStar_luminosity) {
	const double temp = 5772;
	const double solar_lum = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(solar_radius) * m::pow<4>(temp);
	Star sun(temp, solar_radius, 5);
	BOOST_CHECK_CLOSE(sun.luminosity(), solar_lum, 0.1);
}
