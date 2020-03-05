#ifndef FREDDI_ROCHELOBE_HPP
#define FREDDI_ROCHELOBE_HPP

#include <geometry.hpp>


class RochePotential {
public:
	const double mass_ratio;
	const double volume_radius;
public:
	RochePotential(double mass_ratio);
	double omega(double r, double lambda, double nu) const;
	double domega_dr(double r, double lambda, double nu) const;
	double d2omega_dr2(double r, double lambda, double nu) const;
	double d3omega_dr3(double r, double lambda) const;
	double find_r(double omega_value, double lambda, double nu, double guess, double minimum, double maximum) const;
};


class CriticalRocheLobe {
public:
	const RochePotential roche_potential;
	const double mass_ratio;
	const double lagrangian_point;
	const double omega;
	const double polar_radius;
protected:
	static double initializeLagrangianPoint(const RochePotential& roche_potential);
	static double initializePolarRadius(const RochePotential& roche_potential, double lagrangian_point, double omega);
public:
	CriticalRocheLobe(double mass_ratio);
};


class DimensionlessRocheLobe {
public:
	const double mass_ratio;
	const double fill_factor;
	const CriticalRocheLobe critical;
	const double polar_radius;
	const double omega;
public:
	DimensionlessRocheLobe(double mass_ratio, double fill_factor);
	double r(double lambda, double nu) const;
	double r(const Vec3& vec) const;
};


class RocheLobe {
public:
	const double semiaxis;
protected:
	const DimensionlessRocheLobe dimensionless;
public:
	RocheLobe(const double semiaxis, double mass_ratio, double fill_factor);
	double r(const Vec3& vec) const;
};


#endif //FREDDI_ROCHELOBE_HPP
