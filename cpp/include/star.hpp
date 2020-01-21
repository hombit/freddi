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
public:
	IrrSource(const Vec3& position);
	virtual ~IrrSource() = 0;
public:
	const Vec3& position() const;
	double irr_flux(const Vec3& coord, const UnitVec3& normal) const;
	virtual double irr_dLdOmega(const UnitVec3& direction) const = 0;
	virtual bool shadow(const UnitVec3& direction) const = 0; // true for shadowed direction
protected:
	static double cos_object(const UnitVec3& unit_distance, const UnitVec3& normal);
};


class PointLikeSource: virtual public IrrSource {
protected:
	double luminosity_;
public:
	PointLikeSource(const Vec3& position, double luminosity);
	~PointLikeSource() override = 0;
public:
	double luminosity() const;
};


class PointSource: public PointLikeSource {
public:
	using PointLikeSource::PointLikeSource;
	~PointSource() override = 0;
public:
	double irr_dLdOmega(const UnitVec3& direction) const override;
};


class ElementaryPlainSource: public PointLikeSource {
protected:
	UnitVec3 plain_normal_;
public:
	ElementaryPlainSource(const Vec3& position, const UnitVec3& plain_normal, double luminosity);
	~ElementaryPlainSource() override = 0;
public:
	const UnitVec3& plain_normal() const;
	double irr_dLdOmega(const UnitVec3& direction) const override;
};


class DiskShadowSource: virtual public IrrSource {
protected:
	double relative_semiheight_squared_;
public:
	DiskShadowSource(const Vec3& position, double relative_semiheight);
	~DiskShadowSource() override = 0;
public:
	double relative_semiheight_squared() const;
	bool shadow(const UnitVec3& direction) const override;
};


class PointAccretorSource: public PointSource, public DiskShadowSource {
public:
	PointAccretorSource(const Vec3& position, double luminosity, double relative_semiheight);
	~PointAccretorSource() override = default;
};


class CentralDiskSource: public ElementaryPlainSource, public DiskShadowSource {
public:
	CentralDiskSource(const Vec3& position, const UnitVec3& plain_normal,
					  double luminosity, double relative_semiheight);
	~CentralDiskSource() override = default;
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
