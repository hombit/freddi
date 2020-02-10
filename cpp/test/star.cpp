#include <array>
#include <memory>

#include <boost/math/special_functions/pow.hpp>

#include <arguments.hpp>
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
	BOOST_CHECK_CLOSE(sun.luminosity_direction({0.0, 0.0}), solar_lum / (4. * M_PI), 0.1);
}


constexpr static const std::array<std::array<double, 101>, 2> discostar_lum_dir = {{
		{9.14599682e+28, 9.19313081e+28, 9.32588450e+28, 9.55946131e+28,
				9.87902830e+28, 1.03088089e+29, 1.08400609e+29, 1.14741218e+29,
				1.22296603e+29, 1.30816196e+29, 1.40013268e+29, 1.50032581e+29,
				1.60331219e+29, 1.71093048e+29, 1.81876444e+29, 1.92570337e+29,
				2.03361147e+29, 2.15485251e+29, 2.29509555e+29, 2.45414750e+29,
				2.63282268e+29, 2.82916625e+29, 3.03987385e+29, 3.26360110e+29,
				3.50037035e+29, 3.74816288e+29, 4.00306699e+29, 4.26362628e+29,
				4.52781393e+29, 4.79317563e+29, 5.05979607e+29, 5.32681456e+29,
				5.59233080e+29, 5.85449542e+29, 6.11305262e+29, 6.36483210e+29,
				6.60732615e+29, 6.83824325e+29, 7.05610907e+29, 7.26284199e+29,
				7.45421509e+29, 7.63340011e+29, 7.79684288e+29, 7.94547350e+29,
				8.08004678e+29, 8.19714652e+29, 8.29530472e+29, 8.37428786e+29,
				8.43080199e+29, 8.46570242e+29, 8.47708984e+29, 8.46570242e+29,
				8.43080199e+29, 8.37428786e+29, 8.29530472e+29, 8.19714652e+29,
				8.08004678e+29, 7.94547350e+29, 7.79684288e+29, 7.63340011e+29,
				7.45421509e+29, 7.26284199e+29, 7.05610907e+29, 6.83824325e+29,
				6.60732615e+29, 6.36483210e+29, 6.11305262e+29, 5.85449542e+29,
				5.59233080e+29, 5.32681456e+29, 5.05979607e+29, 4.79317563e+29,
				4.52781393e+29, 4.26362628e+29, 4.00306699e+29, 3.74816288e+29,
				3.50037035e+29, 3.26360110e+29, 3.03987385e+29, 2.82916625e+29,
				2.63282268e+29, 2.45414750e+29, 2.29509555e+29, 2.15485251e+29,
				2.03361147e+29, 1.92570337e+29, 1.81876444e+29, 1.71093048e+29,
				1.60331219e+29, 1.50032581e+29, 1.40013268e+29, 1.30816196e+29,
				1.22296603e+29, 1.14741218e+29, 1.08400609e+29, 1.03088089e+29,
				9.87902830e+28, 9.55946131e+28, 9.32588450e+28, 9.19313081e+28,
				9.14599682e+28},
		{1.54082003e+29, 1.55506034e+29, 1.59486509e+29, 1.66591577e+29,
				1.76312151e+29, 1.89581573e+29, 2.06201852e+29, 2.26368509e+29,
				2.50934907e+29, 2.78918841e+29, 3.09074141e+29, 3.42081614e+29,
				3.75490808e+29, 4.10180312e+29, 4.44113766e+29, 4.76851910e+29,
				5.08758125e+29, 5.44191962e+29, 5.86573965e+29, 6.36778033e+29,
				6.95991822e+29, 7.64154941e+29, 8.40096945e+29, 9.23424001e+29,
				1.01436935e+30, 1.11197037e+30, 1.21418172e+30, 1.32032904e+30,
				1.42953504e+30, 1.54073604e+30, 1.65394303e+30, 1.76872650e+30,
				1.88421025e+30, 1.99945069e+30, 2.11421641e+30, 2.22705592e+30,
				2.33664414e+30, 2.44190253e+30, 2.54202987e+30, 2.63784601e+30,
				2.72729397e+30, 2.81169751e+30, 2.88931377e+30, 2.96040365e+30,
				3.02505206e+30, 3.08160646e+30, 3.12917932e+30, 3.16739154e+30,
				3.19485410e+30, 3.21174691e+30, 3.21730283e+30, 3.21174691e+30,
				3.19485410e+30, 3.16739154e+30, 3.12917932e+30, 3.08160646e+30,
				3.02505206e+30, 2.96040365e+30, 2.88931377e+30, 2.81169751e+30,
				2.72729397e+30, 2.63784601e+30, 2.54202987e+30, 2.44190253e+30,
				2.33664414e+30, 2.22705592e+30, 2.11421641e+30, 1.99945069e+30,
				1.88421025e+30, 1.76872650e+30, 1.65394303e+30, 1.54073604e+30,
				1.42953504e+30, 1.32032904e+30, 1.21418172e+30, 1.11197037e+30,
				1.01436935e+30, 9.23424001e+29, 8.40096945e+29, 7.64154941e+29,
				6.95991822e+29, 6.36778033e+29, 5.86573965e+29, 5.44191962e+29,
				5.08758125e+29, 4.76851910e+29, 4.44113766e+29, 4.10180312e+29,
				3.75490808e+29, 3.42081614e+29, 3.09074141e+29, 2.78918841e+29,
				2.50934907e+29, 2.26368509e+29, 2.06201852e+29, 1.89581573e+29,
				1.76312151e+29, 1.66591577e+29, 1.59486509e+29, 1.55506034e+29,
				1.54082003e+29}
}};
constexpr static const std::array<double, 1/*2*/> discostar_lambdas = {{angstromToCm(6410.0)/*, angstromToCm(3465)*/}};

//#include <fstream>
BOOST_AUTO_TEST_CASE(discostar_comparison) {
	const double Mx = sunToGram(1.4);
	const double Mopt = sunToGram(0.55);
	const double Topt = 4000.0;
	const double semiaxis =  3.126121552e11;
	const double Height2R = 0.05;
	const double Lns = 7e36;
	const double Ldisk = 1e37;
	const double albedo = 0.9;
	const double inclination = 40.0 / 180.0 * M_PI;
	const double delta_lambda = angstromToCm(1.0);

	const double Ropt = 0.9 * semiaxis * BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(Mopt / Mx);
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, albedo * Ldisk, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, albedo * Lns, Height2R));
	auto star = IrradiatedStar(std::move(sources), Topt, Ropt, 3);

	const size_t n = 101;
	for (size_t i_lambda = 0; i_lambda < discostar_lambdas.size(); ++i_lambda) {
		const double lambda = discostar_lambdas[i_lambda];
//		std::ofstream file("/Users/hombit/tmp/" + std::to_string(cmToAngstrom(lambda)) + ".dat");
		for (size_t i = 0; i < n; ++i) {
			const double phase = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n - 1);
			const double value = delta_lambda * GSL_CONST_CGSM_SPEED_OF_LIGHT / m::pow<2>(lambda)
								 * star.luminosity_direction({inclination, phase}, lambda);
//			file << phase << "\t" << value << std::endl;
			BOOST_CHECK_CLOSE(value, discostar_lum_dir[i_lambda][i], 15);
		}
	}
}
