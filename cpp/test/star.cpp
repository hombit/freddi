//#define OUTPUT
#define BENCH

#include <array>
#include <memory>
#include <vector>
#ifdef OUTPUT
#include <fstream>
#endif
#ifdef BENCH
#include <chrono>
#include <iostream>
#endif

#include <boost/math/special_functions/pow.hpp>

#include <constants.hpp>
#include <geometry.hpp>
#include <passband.hpp>
#include <rochelobe.hpp>
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


BOOST_AUTO_TEST_CASE(testStar_Roche_vs_spherical_small_roche_fill) {
	const double temp = 5772;
	const double solar_lum =
			GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * 4. * M_PI * m::pow<2>(solar_radius) * m::pow<4>(temp);
	const double radius = solar_radius;
	const double mass_ratio = 1.0;
	const double roche_fill = 0.1;

	const double critical_polar_radius_semiaxis = CriticalRocheLobe(mass_ratio).polar_radius;
	const double semiaxis = radius / critical_polar_radius_semiaxis / roche_fill;
	const RocheLobe roche_lobe(semiaxis, mass_ratio, roche_fill);

	Star sphere(temp, radius, 5);
	Star roche(temp, roche_lobe, 5);

	BOOST_CHECK_CLOSE(sphere.luminosity(), roche.luminosity(), 1);
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
constexpr static const std::array<std::array<double, 101>, 3> discostar_lum_dir_full_roche = {{
	{2.70873836e+28, 2.72218658e+28, 2.76042356e+28, 2.82734835e+28, 2.91911195e+28, 3.04223873e+28, 3.19493535e+28, 3.37807021e+28, 3.59583705e+28, 3.84146081e+28, 4.10896078e+28, 4.39665156e+28, 4.69954712e+28, 5.01734479e+28, 5.34844766e+28, 5.69328473e+28, 6.05208327e+28, 6.42304407e+28, 6.80800743e+28, 7.20437101e+28, 7.61228648e+28, 8.03497607e+28, 8.46800262e+28, 8.90802593e+28, 9.35637493e+28, 9.81051861e+28, 1.02679608e+29, 1.07276000e+29, 1.11883839e+29, 1.16454284e+29, 1.20985611e+29, 1.25487627e+29, 1.29894242e+29, 1.34170908e+29, 1.38311576e+29, 1.42266862e+29, 1.46024988e+29, 1.49555242e+29, 1.52826332e+29, 1.55828778e+29, 1.58541802e+29, 1.60993205e+29, 1.63179608e+29, 1.65142495e+29, 1.66933278e+29, 1.68495723e+29, 1.69808954e+29, 1.70880863e+29, 1.71638353e+29, 1.72114285e+29, 1.72266280e+29, 1.72114285e+29, 1.71638353e+29, 1.70880863e+29, 1.69808954e+29, 1.68495723e+29, 1.66933278e+29, 1.65142495e+29, 1.63179608e+29, 1.60993205e+29, 1.58541802e+29, 1.55828778e+29, 1.52826332e+29, 1.49555242e+29, 1.46024988e+29, 1.42266862e+29, 1.38311576e+29, 1.34170908e+29, 1.29894242e+29, 1.25487627e+29, 1.20985611e+29, 1.16454284e+29, 1.11883839e+29, 1.07276000e+29, 1.02679608e+29, 9.81051861e+28, 9.35637493e+28, 8.90802593e+28, 8.46800262e+28, 8.03497607e+28, 7.61228648e+28, 7.20437101e+28, 6.80800743e+28, 6.42304407e+28, 6.05208327e+28, 5.69328473e+28, 5.34844766e+28, 5.01734479e+28, 4.69954712e+28, 4.39665156e+28, 4.10896078e+28, 3.84146081e+28, 3.59583705e+28, 3.37807021e+28, 3.19493535e+28, 3.04223873e+28, 2.91911195e+28, 2.82734835e+28, 2.76042356e+28, 2.72218658e+28, 2.70873836e+28},
	{3.14419117e+28, 3.16311698e+28, 3.21690693e+28, 3.31127946e+28, 3.44083168e+28, 3.61512825e+28, 3.83184271e+28, 4.09258646e+28, 4.40367538e+28, 4.75535477e+28, 5.13899824e+28, 5.55211479e+28, 5.98712759e+28, 6.44366938e+28, 6.91946251e+28, 7.41527445e+28, 7.93135321e+28, 8.46524126e+28, 9.01984723e+28, 9.59140688e+28, 1.01803447e+29, 1.07916801e+29, 1.14188583e+29, 1.20570394e+29, 1.27083399e+29, 1.33690644e+29, 1.40355406e+29, 1.47062029e+29, 1.53795296e+29, 1.60482282e+29, 1.67120244e+29, 1.73724406e+29, 1.80195268e+29, 1.86480491e+29, 1.92571024e+29, 1.98392674e+29, 2.03926675e+29, 2.09129137e+29, 2.13953078e+29, 2.18384219e+29, 2.22391576e+29, 2.26018798e+29, 2.29261341e+29, 2.32178791e+29, 2.34844905e+29, 2.37175712e+29, 2.39137346e+29, 2.40738452e+29, 2.41871951e+29, 2.42583674e+29, 2.42811505e+29, 2.42583674e+29, 2.41871951e+29, 2.40738452e+29, 2.39137346e+29, 2.37175712e+29, 2.34844905e+29, 2.32178791e+29, 2.29261341e+29, 2.26018798e+29, 2.22391576e+29, 2.18384219e+29, 2.13953078e+29, 2.09129137e+29, 2.03926675e+29, 1.98392674e+29, 1.92571024e+29, 1.86480491e+29, 1.80195268e+29, 1.73724406e+29, 1.67120244e+29, 1.60482282e+29, 1.53795296e+29, 1.47062029e+29, 1.40355406e+29, 1.33690644e+29, 1.27083399e+29, 1.20570394e+29, 1.14188583e+29, 1.07916801e+29, 1.01803447e+29, 9.59140688e+28, 9.01984723e+28, 8.46524126e+28, 7.93135321e+28, 7.41527445e+28, 6.91946251e+28, 6.44366938e+28, 5.98712759e+28, 5.55211479e+28, 5.13899824e+28, 4.75535477e+28, 4.40367538e+28, 4.09258646e+28, 3.83184271e+28, 3.61512825e+28, 3.44083168e+28, 3.31127946e+28, 3.21690693e+28, 3.16311698e+28, 3.14419117e+28},
	{4.06682399e+28, 4.10417077e+28, 4.21019019e+28, 4.39795405e+28, 4.65698289e+28, 5.00932469e+28, 5.45229710e+28, 5.99263439e+28, 6.64681678e+28, 7.39351535e+28, 8.21373106e+28, 9.10217435e+28, 1.00382675e+29, 1.10219820e+29, 1.20473028e+29, 1.31168860e+29, 1.42306856e+29, 1.53839861e+29, 1.65848053e+29, 1.78243114e+29, 1.91052508e+29, 2.04419584e+29, 2.18184708e+29, 2.32237523e+29, 2.46645483e+29, 2.61325916e+29, 2.76191293e+29, 2.91212539e+29, 3.06356183e+29, 3.21439545e+29, 3.36459131e+29, 3.51462545e+29, 3.66196675e+29, 3.80526989e+29, 3.94437860e+29, 4.07748693e+29, 4.20410156e+29, 4.32330529e+29, 4.43395741e+29, 4.53583450e+29, 4.62816792e+29, 4.71230043e+29, 4.78810523e+29, 4.85688838e+29, 4.92025314e+29, 4.97608817e+29, 5.02332487e+29, 5.06190599e+29, 5.08937225e+29, 5.10658910e+29, 5.11213957e+29, 5.10658910e+29, 5.08937225e+29, 5.06190599e+29, 5.02332487e+29, 4.97608817e+29, 4.92025314e+29, 4.85688838e+29, 4.78810523e+29, 4.71230043e+29, 4.62816792e+29, 4.53583450e+29, 4.43395741e+29, 4.32330529e+29, 4.20410156e+29, 4.07748693e+29, 3.94437860e+29, 3.80526989e+29, 3.66196675e+29, 3.51462545e+29, 3.36459131e+29, 3.21439545e+29, 3.06356183e+29, 2.91212539e+29, 2.76191293e+29, 2.61325916e+29, 2.46645483e+29, 2.32237523e+29, 2.18184708e+29, 2.04419584e+29, 1.91052508e+29, 1.78243114e+29, 1.65848053e+29, 1.53839861e+29, 1.42306856e+29, 1.31168860e+29, 1.20473028e+29, 1.10219820e+29, 1.00382675e+29, 9.10217435e+28, 8.21373106e+28, 7.39351535e+28, 6.64681678e+28, 5.99263439e+28, 5.45229710e+28, 5.00932469e+28, 4.65698289e+28, 4.39795405e+28, 4.21019019e+28, 4.10417077e+28, 4.06682399e+28},
}};
constexpr static const std::array<double, 3> discostar_lambdas = {{angstromToCm(6410.0), angstromToCm(5510.0), angstromToCm(3465)}};

void compareWithDiscostar(IrradiatedStar&& star, const std::array<std::array<double, 101>, 3>& discostar_lum_dir) {
#ifdef BENCH
	const auto start = std::chrono::high_resolution_clock::now();
#endif

	const double inclination = 40.0 / 180.0 * M_PI;
	const double delta_lambda = angstromToCm(1.0);

	for (size_t i_lambda = 0; i_lambda < discostar_lambdas.size(); ++i_lambda) {
		const double lambda = discostar_lambdas[i_lambda];
#ifdef OUTPUT
		std::ofstream file("/Users/hombit/tmp/" + std::to_string(cmToAngstrom(lambda)) + ".dat");
#endif
		const size_t n = discostar_lum_dir[i_lambda].size();
		for (size_t i = 0; i < n; ++i) {
			const double phase = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n - 1);
			const UnitVec3 direction(inclination, phase);
			const double value = delta_lambda * GSL_CONST_CGSM_SPEED_OF_LIGHT / m::pow<2>(lambda)
								 * star.luminosity(direction, lambda);
#ifdef OUTPUT
			const double irr_area = star.integrate([&](size_t i_tr) { return static_cast<double>(star.Qirr()[i_tr] > 0.0); }, direction);
			file << phase << "\t" << value << "\t" << irr_area << std::endl;
#endif
			BOOST_CHECK_CLOSE(value, 4.0 * M_PI * discostar_lum_dir[i_lambda][i], 5);
		}
	}

#ifdef OUTPUT
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
#endif

#ifdef BENCH
	const auto end = std::chrono::high_resolution_clock::now();
	std::cout << "BENCH. Light curves computation time for discostar comparison, without star initialisation (N triangles = " << star.triangles().size() << ")"
#ifdef OUTPUT
		<< " (including time of file output)"
#endif // OUTPUT
		<< ": "
		<< std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0 << " ms"
		<< std::endl;
#endif // BENCH
}

BOOST_AUTO_TEST_CASE(testDiscostar_comparison_spherical_geometry_half_roche) {
	const double Mx = sunToGram(1.4);
	const double Mopt = sunToGram(0.55);
	const double Topt = 4000.0;
	const double semiaxis =  3.126121552e11;
	const double polar_radius_semiaxis = 0.2812;
	const double roche_lobe_fill = 0.5;
	const double Height2R = 0.05;
	const double Lns = 7e36;
	const double Ldisk = 1e37;
	const double albedo = 0.9;

	const double Ropt = roche_lobe_fill * polar_radius_semiaxis * semiaxis;
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Ldisk, albedo, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lns, albedo, Height2R));
	IrradiatedStar star(std::move(sources), Topt, Ropt, 3);

	compareWithDiscostar(std::move(star), discostar_lum_dir_half_roche);
}

BOOST_AUTO_TEST_CASE(testDiscostar_comparison_roche_half_roche) {
	const double Mx = sunToGram(1.4);
	const double Mopt = sunToGram(0.55);
	const double Topt = 4000.0;
	const double semiaxis =  3.126121552e11;
	const double roche_lobe_fill = 0.5;
	const double Height2R = 0.05;
	const double Lns = 7e36;
	const double Ldisk = 1e37;
	const double albedo = 0.9;

	const RocheLobe roche_lobe(semiaxis, Mopt / Mx, roche_lobe_fill);
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Ldisk, albedo, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lns, albedo, Height2R));
	IrradiatedStar star(std::move(sources), Topt, roche_lobe, 3);

	compareWithDiscostar(std::move(star), discostar_lum_dir_half_roche);
}

BOOST_AUTO_TEST_CASE(testDiscostar_comparison_full_roche) {
	const double Mx = sunToGram(1.4);
	const double Mopt = sunToGram(0.55);
	const double Topt = 4000.0;
	const double semiaxis =  3.126121552e11;
	const double roche_lobe_fill = 1.0;
	const double Height2R = 0.05;
	const double Lns = 7e36;
	const double Ldisk = 1e37;
	const double albedo = 0.9;

	const RocheLobe roche_lobe(semiaxis, Mopt / Mx, roche_lobe_fill);
	const Vec3 position = {-semiaxis, 0.0, 0.0};
	const UnitVec3 normal(0.0, 0.0);

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, Ldisk, albedo, Height2R));
	sources.push_back(std::make_unique<PointAccretorSource>(position, Lns, albedo, Height2R));
	IrradiatedStar star(std::move(sources), Topt, roche_lobe, 4);

	compareWithDiscostar(std::move(star), discostar_lum_dir_full_roche);
}

#ifdef BENCH
const size_t bench_duration_ms = 3000;

BOOST_AUTO_TEST_CASE(benchRoche_star5) {
	const auto start = std::chrono::high_resolution_clock::now();
	size_t count = 0;

	for (; std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() < bench_duration_ms; ++count) {
		const double temp = 5000;
		const double mass_ratio = 0.1;
		const double roche_fill = 0.8;
		const double semiaxis = 1e12;

		const RocheLobe roche_lobe(semiaxis, mass_ratio, roche_fill);
		Star roche(temp, roche_lobe, 5);
	}

	const auto end = std::chrono::high_resolution_clock::now();

	std::cout
		<< "BENCH. Roche lobe star initialisation time (lod=5): "
		<< std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * count) << " ms"
		<< std::endl;
}
#endif
