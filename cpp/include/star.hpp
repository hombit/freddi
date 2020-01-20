#ifndef FREDDI_STAR_HPP
#define FREDDI_STAR_HPP

#include <functional> // function
#include <optional>
#include <vector>

#include <geometry.hpp>
#include <util.hpp>


class Star {
protected:
	std::vector<Triangle> triangles_;
	vecd Tth_;
	std::optional<vecd> Teff_;
	std::optional<double> luminosity_;
private:
	static std::vector<Triangle> initializeTriangles(double radius, unsigned short grid_scale);
public:
	Star(double temp, double radius, unsigned short grid_scale);
protected:
	double integrate(const vecd& values) const;
	double integrate(std::function<double (size_t)> func) const;
	double integrate(const vecd& values, const UnitVec3& direction) const;
	double integrate(std::function<double (size_t)> func, const UnitVec3& direction) const;
public:
	const std::vector<Triangle>& triangles() const;
	const vecd& Tth() const;
	const vecd& Teff();
	double luminosity();
};


#endif //FREDDI_STAR_HPP
