#include <memory>
#include <vector>

#include <boost/math/special_functions/pow.hpp>

#include <constants.hpp>
#include <geometry.hpp>
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


BOOST_AUTO_TEST_CASE(testStar_orbital_light_curve) {
	const double inclination = M_PI_4;

	const size_t n_phase = 100;
	std::vector<double> phases;
	for (size_t i = 0; i < n_phase; ++i) {
		phases.push_back(2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n_phase));
	}

	std::vector<std::unique_ptr<IrrSource>> sources;
	sources.push_back(std::make_unique<PointAccretorSource>(Vec3(-3.13e11, 0.0, 0.0), 6e37, 0.076));
	sources.push_back(std::make_unique<CentralDiskSource>(Vec3(-2e11, 0.0, 0.0), UnitVec3(0, 0), 4e37, 0.076));
	IrradiatedStar star(std::move(sources), 3000, 9.4e10, 3);
//	std::cout << star.triangles() << std::endl;
//	IrradiatedStar star({}, 4000, 3.5e10, 3);

	for (const auto& phase : phases) {
		const double fluxR = star.dLdOmega({inclination, phase}, angstromToCm(6500)) / m::pow<2>(kpcToCm(5.0));
		const double fluxJ = star.dLdOmega({inclination, phase}, angstromToCm(12400)) / m::pow<2>(kpcToCm(5.0));
		std::cout << phase << "\t" << fluxR << "\t" << fluxJ << std::endl;
	}
}
