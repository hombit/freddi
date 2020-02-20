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
	BOOST_CHECK_CLOSE(sun.luminosity({4.3, 4.8}), solar_lum, 0.1);
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


BOOST_AUTO_TEST_CASE(testIrrStar_point_source) {
	const double temp = 5000;
	const double radius = 5e10;
	const double Lx = 1e36;
	const double semiaxis = 1e12;

	const Vec3 position(-semiaxis, 0.0, 0.0);
	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lx, 0.0, 0.0));
	IrradiatedStar star(std::move(sources), temp, radius, 3);

	const double lum_th = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(radius) * m::pow<4>(temp);
	// https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
	// theta = arcsin(R / a)
	const double omega_star = 2.0 * M_PI * (1.0 - std::sqrt(1.0 - m::pow<2>(radius / semiaxis)));
	const double lum_irr = Lx * omega_star / (4.0 * M_PI);

	BOOST_CHECK_CLOSE(star.luminosity(), lum_th + lum_irr, 1);
}


BOOST_AUTO_TEST_CASE(testIrrStar_albedo) {
	const double radius = 5e10;
	const double semiaxis = 1e12;
	const double Lx = 1e36;

	const Vec3 position(-semiaxis, 0.0, 0.0);
	IrradiatedStar star({}, 0, radius, 3);

	const double albedo0 = 0.0;
	IrradiatedStar::sources_t sources0;
	sources0.push_back(std::make_unique<PointAccretorSource>(position, Lx, albedo0, 0.0));
	star.set_sources(std::move(sources0));
	const double lum0 = star.luminosity();

	const double albedo1 = 0.37;
	IrradiatedStar::sources_t sources1;
	sources1.push_back(std::make_unique<PointAccretorSource>(position, Lx, albedo1, 0.0));
	star.set_sources(std::move(sources1));
	const double lum1 = star.luminosity();

	BOOST_CHECK_CLOSE(lum0 / (1.0 - albedo0), lum1 / (1.0 - albedo1), 1e-6);
}


BOOST_AUTO_TEST_CASE(testIrrStar_shadowed_source) {
	const double temp = 5000;
	const double radius = 5e10;
	const double Lx = 1e100;
	const double semiaxis = 1e12;

	const Vec3 position(-semiaxis, 0.0, 0.0);
	const UnitVec3 normal(0.0, 0.0);
	const double sin = radius / semiaxis;
	// disk h/r = tg = sin / sqrt(1 - sin^2)
	const double height2radius = sin / std::sqrt(1.0 - m::pow<2>(sin));
	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lx, 0.0, height2radius));
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Lx, 0.0, height2radius));
	IrradiatedStar star(std::move(sources), temp, radius, 3);

	const double lum_th = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(radius) * m::pow<4>(temp);

	BOOST_CHECK_CLOSE(star.luminosity(), lum_th, 1);
}


BOOST_AUTO_TEST_CASE(testIrrStar_face_on) {
	const double temp = 5000;
	const double radius = 5e10;
	const double semiaxis = 1e12;
	const double Lns = 1e36;
	const double Ldisk = 1e37;

	const Vec3 position(-semiaxis, 0.0, 0.0);
	const UnitVec3 normal(0.0, 0.0);
	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lns, 0.0, 0.0));
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Ldisk, 0.0, 0.0));
	IrradiatedStar star(std::move(sources), temp, radius, 3);

	const double lum0 = star.luminosity({0.0, 0.0});
	for (size_t phase = 1; phase < 11; ++phase) {
		BOOST_CHECK_EQUAL(star.luminosity({0.0, static_cast<double>(phase)}), lum0);
	}
}


constexpr static const std::array<std::array<double, 101>, 3> discostar_lum_dir_half_roche = {{
	{6.26341288e+27, 6.28321725e+27, 6.34231379e+27, 6.43962943e+27, 6.57514164e+27, 6.74818369e+27, 6.95734909e+27, 7.20231594e+27, 7.48252085e+27, 7.79682031e+27, 8.14400425e+27, 8.52332331e+27, 8.93418145e+27, 9.37613526e+27, 9.84866355e+27, 1.03511280e+28, 1.08832552e+28, 1.14443430e+28, 1.20333919e+28, 1.26489893e+28, 1.32896451e+28, 1.39541849e+28, 1.46409376e+28, 1.53473327e+28, 1.60709322e+28, 1.68095047e+28, 1.75603341e+28, 1.83202585e+28, 1.90856546e+28, 1.98531519e+28, 2.06195814e+28, 2.13809670e+28, 2.21329500e+28, 2.28716231e+28, 2.35931589e+28, 2.42934953e+28, 2.49684702e+28, 2.56137104e+28, 2.62254427e+28, 2.68000300e+28, 2.73340265e+28, 2.78244240e+28, 2.82686660e+28, 2.86649693e+28, 2.90118073e+28, 2.93074079e+28, 2.95501008e+28, 2.97396520e+28, 2.98754452e+28, 2.99565829e+28, 2.99835849e+28, 2.99565829e+28, 2.98754452e+28, 2.97396520e+28, 2.95501008e+28, 2.93074079e+28, 2.90118073e+28, 2.86649693e+28, 2.82686660e+28, 2.78244240e+28, 2.73340265e+28, 2.68000300e+28, 2.62254427e+28, 2.56137104e+28, 2.49684702e+28, 2.42934953e+28, 2.35931589e+28, 2.28716231e+28, 2.21329500e+28, 2.13809670e+28, 2.06195814e+28, 1.98531519e+28, 1.90856546e+28, 1.83202585e+28, 1.75603341e+28, 1.68095047e+28, 1.60709322e+28, 1.53473327e+28, 1.46409376e+28, 1.39541849e+28, 1.32896451e+28, 1.26489893e+28, 1.20333919e+28, 1.14443430e+28, 1.08832552e+28, 1.03511280e+28, 9.84866355e+27, 9.37613526e+27, 8.93418145e+27, 8.52332331e+27, 8.14400425e+27, 7.79682031e+27, 7.48252085e+27, 7.20231594e+27, 6.95734909e+27, 6.74818369e+27, 6.57514164e+27, 6.43962943e+27, 6.34231379e+27, 6.28321725e+27, 6.26341288e+27},
	{7.19060564e+27, 7.21820361e+27, 7.30056324e+27, 7.43627927e+27, 7.62542539e+27, 7.86718454e+27, 8.15971625e+27, 8.50271408e+27, 8.89560639e+27, 9.33696592e+27, 9.82522408e+27, 1.03594499e+28, 1.09389815e+28, 1.15631925e+28, 1.22312570e+28, 1.29423266e+28, 1.36960860e+28, 1.44915528e+28, 1.53273268e+28, 1.62014953e+28, 1.71119808e+28, 1.80571200e+28, 1.90345882e+28, 2.00407808e+28, 2.10722410e+28, 2.21257850e+28, 2.31975467e+28, 2.42830045e+28, 2.53769878e+28, 2.64746650e+28, 2.75714509e+28, 2.86616118e+28, 2.97389153e+28, 3.07977374e+28, 3.18325136e+28, 3.28374210e+28, 3.38065013e+28, 3.47334392e+28, 3.56127882e+28, 3.64394402e+28, 3.72084723e+28, 3.79154029e+28, 3.85564460e+28, 3.91288986e+28, 3.96303781e+28, 4.00581235e+28, 4.04095984e+28, 4.06843022e+28, 4.08812214e+28, 4.09989694e+28, 4.10381619e+28, 4.09989694e+28, 4.08812214e+28, 4.06843022e+28, 4.04095984e+28, 4.00581235e+28, 3.96303781e+28, 3.91288986e+28, 3.85564460e+28, 3.79154029e+28, 3.72084723e+28, 3.64394402e+28, 3.56127882e+28, 3.47334392e+28, 3.38065013e+28, 3.28374210e+28, 3.18325136e+28, 3.07977374e+28, 2.97389153e+28, 2.86616118e+28, 2.75714509e+28, 2.64746650e+28, 2.53769878e+28, 2.42830045e+28, 2.31975467e+28, 2.21257850e+28, 2.10722410e+28, 2.00407808e+28, 1.90345882e+28, 1.80571200e+28, 1.71119808e+28, 1.62014953e+28, 1.53273268e+28, 1.44915528e+28, 1.36960860e+28, 1.29423266e+28, 1.22312570e+28, 1.15631925e+28, 1.09389815e+28, 1.03594499e+28, 9.82522408e+27, 9.33696592e+27, 8.89560639e+27, 8.50271408e+27, 8.15971625e+27, 7.86718454e+27, 7.62542539e+27, 7.43627927e+27, 7.30056324e+27, 7.21820361e+27, 7.19060564e+27},
	{8.25610088e+27, 8.30556498e+27, 8.45327410e+27, 8.69743202e+27, 9.03936038e+27, 9.47874405e+27, 1.00134052e+28, 1.06442695e+28, 1.13721689e+28, 1.21960216e+28, 1.31140030e+28, 1.41255971e+28, 1.52309752e+28, 1.64293134e+28, 1.77184761e+28, 1.90970485e+28, 2.05648536e+28, 2.21200734e+28, 2.37599321e+28, 2.54810494e+28, 2.72794491e+28, 2.91519085e+28, 3.10940349e+28, 3.30989182e+28, 3.51596929e+28, 3.72699873e+28, 3.94220665e+28, 4.16068437e+28, 4.38139354e+28, 4.60335844e+28, 4.82562762e+28, 5.04702417e+28, 5.26628686e+28, 5.48226875e+28, 5.69380830e+28, 5.89973447e+28, 6.09883708e+28, 6.28980401e+28, 6.47151432e+28, 6.64299149e+28, 6.80320493e+28, 6.95110600e+28, 7.08580032e+28, 7.20661362e+28, 7.31289003e+28, 7.40387795e+28, 7.47891472e+28, 7.53774858e+28, 7.58004392e+28, 7.60540529e+28, 7.61385393e+28, 7.60540529e+28, 7.58004392e+28, 7.53774858e+28, 7.47891472e+28, 7.40387795e+28, 7.31289003e+28, 7.20661362e+28, 7.08580032e+28, 6.95110600e+28, 6.80320493e+28, 6.64299149e+28, 6.47151432e+28, 6.28980401e+28, 6.09883708e+28, 5.89973447e+28, 5.69380830e+28, 5.48226875e+28, 5.26628686e+28, 5.04702417e+28, 4.82562762e+28, 4.60335844e+28, 4.38139354e+28, 4.16068437e+28, 3.94220665e+28, 3.72699873e+28, 3.51596929e+28, 3.30989182e+28, 3.10940349e+28, 2.91519085e+28, 2.72794491e+28, 2.54810494e+28, 2.37599321e+28, 2.21200734e+28, 2.05648536e+28, 1.90970485e+28, 1.77184761e+28, 1.64293134e+28, 1.52309752e+28, 1.41255971e+28, 1.31140030e+28, 1.21960216e+28, 1.13721689e+28, 1.06442695e+28, 1.00134052e+28, 9.47874405e+27, 9.03936038e+27, 8.69743202e+27, 8.45327410e+27, 8.30556498e+27, 8.25610088e+27},
}};
constexpr static const std::array<double, 3> discostar_lambdas = {{angstromToCm(6410.0), angstromToCm(5510.0), angstromToCm(3465)}};

//#include <fsteam>
BOOST_AUTO_TEST_CASE(discostar_comparison) {
	const double Mx = sunToGram(1.4);
	const double Mopt = sunToGram(0.55);
	const double Topt = 4000.0;
	const double semiaxis =  3.126121552e11;
	const double roche_fill = 0.5;
	const double polar_radius_semiaxis = 0.2812;
	const double Height2R = 0.05;
	const double Lns = 7e36;
	const double Ldisk = 1e37;
	const double albedo = 0.9;
	const double inclination = 40.0 / 180.0 * M_PI;
	const double delta_lambda = angstromToCm(1.0);

//	const double Ropt = semiaxis * BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(Mopt / Mx);
	const double Ropt = roche_fill * polar_radius_semiaxis * semiaxis;
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Ldisk, albedo, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lns, albedo, Height2R));
	IrradiatedStar star(std::move(sources), Topt, Ropt, 3);

	for (size_t i_lambda = 0; i_lambda < discostar_lambdas.size(); ++i_lambda) {
		const double lambda = discostar_lambdas[i_lambda];
//		std::ofstream file("/Users/hombit/tmp/" + std::to_string(cmToAngstrom(lambda)) + ".dat");
		const size_t n = discostar_lum_dir_half_roche[i_lambda].size();
		for (size_t i = 0; i < n; ++i) {
			const double phase = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n - 1);
			const UnitVec3 direction(inclination, phase);
			const double value = delta_lambda * GSL_CONST_CGSM_SPEED_OF_LIGHT / m::pow<2>(lambda)
								 * star.luminosity(direction, lambda);
//			const double irr_area = star.integrate([&](size_t i_tr) { return static_cast<double>(star.Qirr()[i_tr] > 0.0); }, direction);
//			file << phase << "\t" << value << "\t" << irr_area << std::endl;
			BOOST_CHECK_CLOSE(value, 4.0 * M_PI * discostar_lum_dir_half_roche[i_lambda][i], 5);
		}
	}

	/*
	std::ofstream file("/Users/hombit/tmp/star.dat");
	file << "vertex0_x\tvertex0_y\tvertex0_z\tvertex1_x\tvertex1_y\tvertex1_z\tvertex2_x\tvertex2_y\tvertex2_z\tcenter_x\tcenter_y\tcenter_z\tTeff\n";
	for (size_t i_tr = 0; i_tr < star.triangles().size(); ++i_tr) {
		const auto& tr = star.triangles()[i_tr];
		for (const auto& vertex : tr.vertices()) {
			file << vertex.x() << "\t" << vertex.y() << "\t" << vertex.z() << "\t";
		}
		const auto& center = star.triangles()[i_tr].center();
		file << center.x() << "\t" << center.y() << "\t" << center.z() << "\t";
		file << star.Teff()[i_tr] << std::endl;
	}
	*/
}
