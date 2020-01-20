#ifndef FREDDI_STAR_HPP
#define FREDDI_STAR_HPP

#include <functional> // function
#include <memory> // unique_ptr
#include <optional>
#include <vector>

#include <arguments.hpp>
#include <geometry.hpp>
#include <passband.hpp>
#include <util.hpp>


class Star {
protected:
	struct IrradiatedProperties {
		std::optional<vald> Qirr;
		std::optional<vald> Teff;
		std::optional<double> luminosity;
	};
protected:
	const std::vector<Triangle> triangles_;
	const vald Tth_;
private:
	static std::vector<Triangle> initializeTriangles(double radius, unsigned short grid_scale);
protected:
	IrradiatedProperties irr_;
protected:
	void invalidate_irradiated_properties();
public:
	Star(double temp, double radius, unsigned short lod);
protected:
	double integrate(std::function<double (size_t)>&& func) const;
	double integrate(std::function<double (size_t)>&& func, const UnitVec3& direction) const;
public:
	const std::vector<Triangle>& triangles() const;
	const vald& Tth() const;
	virtual const vald& Qirr();
	const vald& Teff();
	double luminosity();
	double dLdOmega(const UnitVec3& direction);
	double dLdOmega(const UnitVec3& direction, double lambda);
	double dLdOmega(const UnitVec3& direction, const Passband& passband);
};


class IrrSource {
protected:
	Vec3 position_;
	double luminosity_;
public:
	IrrSource(const Vec3& position);
	virtual double irr_flux(const Vec3& coord, const UnitVec3& normal) const;
	const Vec3& position() const;
	double luminosity() const;
	void set_luminosity(double value);
	static double cos_object(const UnitVec3& unit_distance, const UnitVec3& normal);
};


class PointSource: public IrrSource {
public:
	using IrrSource::IrrSource;
public:
	double irr_flux(const Vec3& coord, const UnitVec3& normal) const override;
};


class PlainSource: public IrrSource {
protected:
	UnitVec3 plain_normal_;
public:
	PlainSource(const Vec3& position, const UnitVec3& plain_normal);
public:
	double irr_flux(const Vec3& coord, const UnitVec3& normal) const override;
	const UnitVec3& plain_normal() const;
};


class IrradiatedStar: public Star {
protected:
	std::vector<std::unique_ptr<IrrSource>> sources_;
public:
	IrradiatedStar(std::vector<std::unique_ptr<IrrSource>>&& sources,
			double temp, double radius, unsigned short lod);
public:
	const std::vector<std::unique_ptr<IrrSource>>& sources() const;
	const vald& Qirr() override;
};


#endif //FREDDI_STAR_HPP
