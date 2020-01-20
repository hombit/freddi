#ifndef FREDDI_STAR_HPP
#define FREDDI_STAR_HPP

#include <functional> // function
#include <optional>
#include <vector>

#include <arguments.hpp>
#include <geometry.hpp>
#include <util.hpp>


class Star {
protected:
	struct IrradiatedProperties {
		std::optional<vald> Tirr;
		std::optional<vald> Teff;
		std::optional<double> luminosity;
	};
protected:
	const std::vector<Triangle> triangles_;
	const vald Tth_;
private:
	static std::vector<Triangle> initializeTriangles(double radius, unsigned short grid_scale);
	IrradiatedProperties irr_;
protected:
	void invalidate_irradiated_properties();
public:
	Star(double temp, double radius, unsigned short grid_scale);
protected:
	double integrate(std::function<double (size_t)> func) const;
	double integrate(std::function<double (size_t)> func, const UnitVec3& direction) const;
public:
	const std::vector<Triangle>& triangles() const;
	const vald& Tth() const;
	virtual const vald& Tirr();
	const vald& Teff();
	double luminosity();
};


class FreddiStar: public Star {
protected:
	FreddiArguments args;
public:
	FreddiStar(const FreddiArguments& args);
};


#endif //FREDDI_STAR_HPP
