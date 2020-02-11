#include <array>
#include <memory>
#include <vector>

#include <boost/math/special_functions/pow.hpp>

#include <arguments.hpp>
#include <constants.hpp>
#include <geometry.hpp>
#include <passband.hpp>
#include <star.hpp>
#include <unit_transformation.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testStar

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>


namespace m = boost::math;


BOOST_AUTO_TEST_CASE(testStar_luminosity) {
	const double temp = 5772;
	const double solar_lum =
			GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(solar_radius) * m::pow<4>(temp);
	Star sun(temp, solar_radius, 5);
	BOOST_CHECK_CLOSE(sun.luminosity(), solar_lum, 0.1);
	BOOST_CHECK_CLOSE(sun.luminosity({0.0, 0.0}), solar_lum, 0.1);
}


BOOST_AUTO_TEST_CASE(testStar_luminosity_passband) {
	const double temp = 5772;
	const double solar_lum =
			GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(solar_radius) * m::pow<4>(temp);
	Star sun(temp, solar_radius, 5);

	std::vector<Passband::PassbandPoint> passband_data;
	const double lambda_min = angstromToCm(100);
	const double lambda_max = angstromToCm(100000);
	const double lambda_count = 1000;
	for (size_t i = 0; i < lambda_count; ++i) {
		const double lambda = lambda_min * std::pow(lambda_max / lambda_min, static_cast<double>(i) / static_cast<double>(lambda_count - 1));
		passband_data.emplace_back(lambda, 1.0);
	}
	const Passband white_filter("white", passband_data);
	BOOST_CHECK_CLOSE(sun.luminosity({0.0, 0.0}, white_filter) * white_filter.t_dnu, solar_lum, 1);

	const double lambda_V = angstromToCm(5510);
	const Passband narrow_filter("V", {{lambda_V, 1.0}, {lambda_V * 1.0001, 1.0}});
	BOOST_CHECK_CLOSE(sun.luminosity({0.0, 0.0}, narrow_filter), sun.luminosity({0.0, 0.0}, lambda_V), 0.01);
}


BOOST_AUTO_TEST_CASE(testIrrPoint_Star_luminosity) {
	const double temp = 5000;
	const double radius = 5e10;
	const double Lx = 1e36;
	const double semiaxis = 1e12;

	const Vec3 position(-semiaxis, 0.0, 0.0);
	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lx, 0.0));
	IrradiatedStar star(std::move(sources), temp, radius, 3);

	const double lum_cold = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(radius) * m::pow<4>(temp);
	// https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
	const double omega_star = 4.0 * M_PI * m::pow<2>(radius / (2.0 * semiaxis));
	const double lum_irr = Lx * omega_star / (4.0 * M_PI);

	BOOST_CHECK_CLOSE(star.luminosity(), lum_cold + lum_irr, 1.0);
}


constexpr static const std::array<std::array<double, 101>, 3> discostar_lum_dir = {{
	{9.14599682e+28, 9.19313081e+28, 9.32586050e+28, 9.55873983e+28, 9.87503788e+28, 1.02988388e+29, 1.08206860e+29, 1.14405555e+29, 1.21716992e+29, 1.29815914e+29, 1.38401410e+29, 1.47354821e+29, 1.56420494e+29, 1.65545214e+29, 1.74618100e+29, 1.83601889e+29, 1.92680932e+29, 2.02989108e+29, 2.15060605e+29, 2.28976863e+29, 2.44767698e+29, 2.62235734e+29, 2.81016323e+29, 3.00986533e+29, 3.22083918e+29, 3.44047759e+29, 3.66685637e+29, 3.89858064e+29, 4.13466441e+29, 4.37287406e+29, 4.61266928e+29, 4.85385999e+29, 5.09346529e+29, 5.32952894e+29, 5.56110272e+29, 5.78565935e+29, 6.00202967e+29, 6.20825141e+29, 6.40241293e+29, 6.58341550e+29, 6.74968027e+29, 6.90176825e+29, 7.03900365e+29, 7.16268277e+29, 7.27451003e+29, 7.37179018e+29, 7.45334566e+29, 7.51914249e+29, 7.56600138e+29, 7.59504002e+29, 7.60446784e+29, 7.59504002e+29, 7.56600138e+29, 7.51914249e+29, 7.45334566e+29, 7.37179018e+29, 7.27451003e+29, 7.16268277e+29, 7.03900365e+29, 6.90176825e+29, 6.74968027e+29, 6.58341550e+29, 6.40241293e+29, 6.20825141e+29, 6.00202967e+29, 5.78565935e+29, 5.56110272e+29, 5.32952894e+29, 5.09346529e+29, 4.85385999e+29, 4.61266928e+29, 4.37287406e+29, 4.13466441e+29, 3.89858064e+29, 3.66685637e+29, 3.44047759e+29, 3.22083918e+29, 3.00986533e+29, 2.81016323e+29, 2.62235734e+29, 2.44767698e+29, 2.28976863e+29, 2.15060605e+29, 2.02989108e+29, 1.92680932e+29, 1.83601889e+29, 1.74618100e+29, 1.65545214e+29, 1.56420494e+29, 1.47354821e+29, 1.38401410e+29, 1.29815914e+29, 1.21716992e+29, 1.14405555e+29, 1.08206860e+29, 1.02988388e+29, 9.87503788e+28, 9.55873983e+28, 9.32586050e+28, 9.19313081e+28, 9.14599682e+28},
	{1.08569751e+29, 1.09242620e+29, 1.11134625e+29, 1.14464596e+29, 1.18990240e+29, 1.25072595e+29, 1.32580281e+29, 1.41523934e+29, 1.52108059e+29, 1.63844079e+29, 1.76277486e+29, 1.89223485e+29, 2.02279059e+29, 2.15365447e+29, 2.28314316e+29, 2.41072062e+29, 2.53889372e+29, 2.68466275e+29, 2.85709328e+29, 3.05808652e+29, 3.28868857e+29, 3.54612838e+29, 3.82484545e+29, 4.12298803e+29, 4.43947762e+29, 4.77021033e+29, 5.11221661e+29, 5.46334676e+29, 5.82209279e+29, 6.18503297e+29, 6.55130191e+29, 6.92057450e+29, 7.28823699e+29, 7.65120301e+29, 8.00793575e+29, 8.35448738e+29, 8.68893051e+29, 9.00823402e+29, 9.30936849e+29, 9.59056231e+29, 9.84929159e+29, 1.00864496e+30, 1.03009514e+30, 1.04946404e+30, 1.06699572e+30, 1.08226902e+30, 1.09508654e+30, 1.10542448e+30, 1.11279959e+30, 1.11736609e+30, 1.11885199e+30, 1.11736609e+30, 1.11279959e+30, 1.10542448e+30, 1.09508654e+30, 1.08226902e+30, 1.06699572e+30, 1.04946404e+30, 1.03009514e+30, 1.00864496e+30, 9.84929159e+29, 9.59056231e+29, 9.30936849e+29, 9.00823402e+29, 8.68893051e+29, 8.35448738e+29, 8.00793575e+29, 7.65120301e+29, 7.28823699e+29, 6.92057450e+29, 6.55130191e+29, 6.18503297e+29, 5.82209279e+29, 5.46334676e+29, 5.11221661e+29, 4.77021033e+29, 4.43947762e+29, 4.12298803e+29, 3.82484545e+29, 3.54612838e+29, 3.28868857e+29, 3.05808652e+29, 2.85709328e+29, 2.68466275e+29, 2.53889372e+29, 2.41072062e+29, 2.28314316e+29, 2.15365447e+29, 2.02279059e+29, 1.89223485e+29, 1.76277486e+29, 1.63844079e+29, 1.52108059e+29, 1.41523934e+29, 1.32580281e+29, 1.25072595e+29, 1.18990240e+29, 1.14464596e+29, 1.11134625e+29, 1.09242620e+29, 1.08569751e+29},
	{1.54082003e+29, 1.55506034e+29, 1.59486402e+29, 1.66584111e+29, 1.76259452e+29, 1.89439675e+29, 2.05898172e+29, 2.25780681e+29, 2.49699888e+29, 2.76369667e+29, 3.04567608e+29, 3.33770254e+29, 3.62701979e+29, 3.91182161e+29, 4.18708210e+29, 4.45159823e+29, 4.70889846e+29, 4.99769441e+29, 5.35099246e+29, 5.78188630e+29, 6.30113673e+29, 6.90609676e+29, 7.58268071e+29, 8.32621467e+29, 9.13213075e+29, 9.98761178e+29, 1.08842135e+30, 1.18160714e+30, 1.27791699e+30, 1.37640262e+30, 1.47680669e+30, 1.57901710e+30, 1.68170554e+30, 1.78391480e+30, 1.88513822e+30, 1.98420804e+30, 2.08045312e+30, 2.17297089e+30, 2.26080464e+30, 2.34340016e+30, 2.41991385e+30, 2.49064164e+30, 2.55518315e+30, 2.61390344e+30, 2.66730765e+30, 2.71408731e+30, 2.75349410e+30, 2.78525201e+30, 2.80803689e+30, 2.82210168e+30, 2.82671204e+30, 2.82210168e+30, 2.80803689e+30, 2.78525201e+30, 2.75349410e+30, 2.71408731e+30, 2.66730765e+30, 2.61390344e+30, 2.55518315e+30, 2.49064164e+30, 2.41991385e+30, 2.34340016e+30, 2.26080464e+30, 2.17297089e+30, 2.08045312e+30, 1.98420804e+30, 1.88513822e+30, 1.78391480e+30, 1.68170554e+30, 1.57901710e+30, 1.47680669e+30, 1.37640262e+30, 1.27791699e+30, 1.18160714e+30, 1.08842135e+30, 9.98761178e+29, 9.13213075e+29, 8.32621467e+29, 7.58268071e+29, 6.90609676e+29, 6.30113673e+29, 5.78188630e+29, 5.35099246e+29, 4.99769441e+29, 4.70889846e+29, 4.45159823e+29, 4.18708210e+29, 3.91182161e+29, 3.62701979e+29, 3.33770254e+29, 3.04567608e+29, 2.76369667e+29, 2.49699888e+29, 2.25780681e+29, 2.05898172e+29, 1.89439675e+29, 1.76259452e+29, 1.66584111e+29, 1.59486402e+29, 1.55506034e+29, 1.54082003e+29}
}};
constexpr static const std::array<double, 3> discostar_lambdas = {{angstromToCm(6410.0), angstromToCm(5510.0), angstromToCm(3465)}};

#include <fstream>
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

	const double Ropt = semiaxis * BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(Mopt / Mx);
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, (1.0 - albedo) * Ldisk, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, (1.0 - albedo) * Lns, Height2R));
	IrradiatedStar star(std::move(sources), Topt, Ropt, 3);

	for (size_t i_lambda = 0; i_lambda < discostar_lambdas.size(); ++i_lambda) {
		const double lambda = discostar_lambdas[i_lambda];
		std::ofstream file("/Users/hombit/tmp/" + std::to_string(cmToAngstrom(lambda)) + ".dat");
		const size_t n = discostar_lum_dir[i_lambda].size();
		for (size_t i = 0; i < n; ++i) {
			const double phase = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n - 1);
			const double value = delta_lambda * GSL_CONST_CGSM_SPEED_OF_LIGHT / m::pow<2>(lambda)
								 * star.luminosity({inclination, phase}, lambda);
			file << phase << "\t" << value << std::endl;
			BOOST_CHECK_CLOSE(value, 4.0 * M_PI * discostar_lum_dir[i_lambda][i], 15);
		}
	}
}
