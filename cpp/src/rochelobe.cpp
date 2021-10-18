#include <algorithm> // min

#include <boost/math/tools/roots.hpp>

#include <orbit.hpp>
#include <rochelobe.hpp>
#include <util.hpp>

using boost::math::tools::halley_iterate;


RochePotential::RochePotential(const double mass_ratio):
		mass_ratio(mass_ratio),
		volume_radius(BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(mass_ratio)) {}

double RochePotential::omega(const double r, const double lambda, const double nu) const {
	// Cherepashchuk, Fismatlit, Moscow, 2013, Volume 1, Eq. 365, but our mass_ratio = 1/q
	return mass_ratio / r
		+ (1.0 / std::sqrt(1.0 + m::pow<2>(r) - 2.0 * r * lambda) - r * lambda)
		+ 0.5 * (1.0 + mass_ratio) * m::pow<2>(r) * (1.0 - m::pow<2>(nu));
}

double RochePotential::domega_dr(const double r, const double lambda, const double nu) const {
	return -mass_ratio / m::pow<2>(r)
		   - ((r - lambda) * std::pow(1.0 + m::pow<2>(r) - 2.0 * r * lambda, -1.5) + lambda)
		   + (1.0 + mass_ratio) * r * (1.0 - m::pow<2>(nu));
}

double RochePotential::d2omega_dr2(const double r, const double lambda, const double nu) const {
	const double rcm2 = 1.0 / std::sqrt(1.0 + m::pow<2>(r) - 2.0 * r * lambda);
	return 2.0 * mass_ratio / m::pow<3>(r)
		   - m::pow<3>(rcm2)
		   + 3.0 * m::pow<2>(r - lambda) * m::pow<5>(rcm2)
		   + (1.0 + mass_ratio) * (1.0 - m::pow<2>(nu));
}

double RochePotential::d3omega_dr3(const double r, const double lambda) const {
	const double rcm2 = 1.0 / std::sqrt(1.0 + m::pow<2>(r) - 2.0 * r * lambda);
	return -6.0 * mass_ratio / m::pow<4>(r)
		+ 9.0 * (r - lambda) * m::pow<5>(rcm2)
		- 15.0 * m::pow<3>(r - lambda) * m::pow<7>(rcm2);
}

double RochePotential::find_r(const double omega_value, const double lambda, const double nu,
							  const double guess, const double minimum, const double maximum) const {
	return halley_iterate(
			[this, omega_value, lambda, nu](double x) {
				return std::make_tuple(omega(x, lambda, nu) - omega_value,
									   domega_dr(x, lambda, nu),
									   d2omega_dr2(x, lambda, nu));
			},
			guess, minimum, maximum,
			std::numeric_limits<double>::digits / 2
	);
}


CriticalRocheLobe::CriticalRocheLobe(const double mass_ratio):
		roche_potential(mass_ratio),
		mass_ratio(mass_ratio),
		lagrangian_point(initializeLagrangianPoint(roche_potential)),
		omega(roche_potential.omega(lagrangian_point, 1.0, 0.0)),
		polar_radius(initializePolarRadius(roche_potential, lagrangian_point, omega)) {}

double CriticalRocheLobe::initializeLagrangianPoint(const RochePotential& roche_potential) {
	// q <= 1
	double q = roche_potential.mass_ratio <= 1.0 ? roche_potential.mass_ratio : 1.0 / roche_potential.mass_ratio;
	// Initial guess is a Hill radius for small mass_ratio, one minus Hill radius for large mass_ratio and 0.5 for
	// mass_ratio = 1
	double x0 = std::cbrt(q / 3.0) + (0.5 - std::cbrt(1.0/3.0)) * q;
	double x1 = 0.5 * x0;
	double x2 = std::min(0.5, 2 * x0);
	if (roche_potential.mass_ratio > 1.0) {
		x0 = 1.0 - x0;
		const double tmp = 1.0 - x1;
		x1 = 1.0 - x2;
		x2 = tmp;
	}

	const double lagrangian_point = halley_iterate(
			[&roche_potential](double x) {
				return std::make_tuple(roche_potential.domega_dr(x, 1.0, 0.0),
						roche_potential.d2omega_dr2(x, 1.0, 0.0),
						roche_potential.d3omega_dr3(x, 1.0));
			},
			x0, x1, x2,
			std::numeric_limits<double>::digits / 2
	);
	return lagrangian_point;
}

double CriticalRocheLobe::initializePolarRadius(const RochePotential& roche_potential, const double lagrangian_point, const double omega) {
	const double z0 = roche_potential.volume_radius;
	const double zmin = 0.5 * z0;
	const double zmax = lagrangian_point;
	return roche_potential.find_r(omega, 0.0, 1.0, z0, zmin, zmax);
}


DimensionlessRocheLobe::DimensionlessRocheLobe(const double mass_ratio, const double fill_factor):
		mass_ratio(mass_ratio),
		fill_factor(fill_factor),
		critical(mass_ratio),
		polar_radius(fill_factor * critical.polar_radius),
		omega(critical.roche_potential.omega(polar_radius, 0.0, 1.0)) {}

double DimensionlessRocheLobe::r(const double lambda, const double nu) const {
	const double r0 = polar_radius;
	const double rmin = 0.5 * r0;
	const double rmax = std::min(2.0 * r0, critical.lagrangian_point);
	return critical.roche_potential.find_r(omega, lambda, nu, r0, rmin, rmax);
}

double DimensionlessRocheLobe::r(const Vec3& vec) const {
	const double lambda = -vec.x() / vec.length();
	const double nu = vec.z() / vec.length();
	return r(lambda, nu);
}


RocheLobe::RocheLobe(const double semiaxis, double mass_ratio, double fill_factor):
		semiaxis(semiaxis),
		dimensionless(mass_ratio, fill_factor) {}

double RocheLobe::r(const Vec3& vec) const {
	return semiaxis * dimensionless.r(vec);
}
