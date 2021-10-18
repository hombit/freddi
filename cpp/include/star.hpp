#ifndef FREDDI_STAR_HPP
#define FREDDI_STAR_HPP

#include <functional> // function
#include <memory> // unique_ptr
#include <optional>
#include <vector>

#include <arguments.hpp>
#include <geometry.hpp>
#include <passband.hpp>
#include <rochelobe.hpp>
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
	static std::vector<Triangle> initializeSphereTriangles(double radius, unsigned short grid_scale);
	static std::vector<Triangle> initializeRocheTriangles(const RocheLobe& roche_lobe, unsigned short grid_scale);
protected:
	IrradiatedProperties irr_;
protected:
	void invalidate_irradiated_properties();
public:
	Star(double temp, double radius, unsigned short lod);
	Star(double temp, const RocheLobe& roche_lobe, unsigned short lod);
	const std::vector<Triangle>& triangles() const;
	const vald& Tth() const;
	virtual const vald& Qirr();
	const vald& Teff();
	double integrate(std::function<double (size_t)>&& func) const;
	double integrate(std::function<double (size_t)>&& func, const UnitVec3& direction) const;
	double luminosity();
	double luminosity(const UnitVec3& direction); // erg/s
	// "Luminosity in direction" is luminosity of isotropic spherical source of radius R with luminosity given by
	// 4\pi R^2 \pi B_\nu. "Luminosity in direction" = 4 \pi \int \dS cos{n} B_nu. Observable flux is
	// "luminosity in direction" / (4 \pi d^2)
	double luminosity(const UnitVec3& direction, double lambda); // erg/s/Hz
	double luminosity(const UnitVec3& direction, const Passband& passband); // erg/s/Hz
};


class IrrSource {
public:
	virtual ~IrrSource() = 0;
public:
	virtual IrrSource* clone() const = 0;
	virtual double irr_flux(const Vec3& coord, const UnitVec3& normal) const = 0;
	virtual bool shadow(const UnitVec3& direction) const = 0; // true for shadowed direction
	virtual double albedo(double cos_object) const = 0; // A part of irradiation reflected by an object
protected:
	static double cos_object(const UnitVec3& direction, const UnitVec3& normal);
};


class PointLikeSource: virtual public IrrSource {
protected:
	Vec3 position_;
	double luminosity_;
public:
	PointLikeSource(const Vec3& position, double luminosity);
	~PointLikeSource() override = 0;
	double irr_flux(const Vec3& coord, const UnitVec3& normal) const override;
	virtual double irr_luminosity(const UnitVec3& direction) const = 0; // "Luminosity in direction"
public:
	const Vec3& position() const;
	double luminosity() const;
};


class PointSource: public PointLikeSource {
public:
	using PointLikeSource::PointLikeSource;
	~PointSource() override = 0;
public:
	double irr_luminosity(const UnitVec3& direction) const override;
};


class ElementaryPlainSource: public PointLikeSource {
protected:
	UnitVec3 plain_normal_;
public:
	ElementaryPlainSource(const Vec3& position, const UnitVec3& plain_normal, double luminosity);
	~ElementaryPlainSource() override = 0;
public:
	const UnitVec3& plain_normal() const;
	double irr_luminosity(const UnitVec3& direction) const override;
};


class ConstantAlbedoSource: virtual public IrrSource {
protected:
	double albedo_;
public:
	ConstantAlbedoSource(double albedo);
	virtual ~ConstantAlbedoSource() = 0;
public:
	double albedo(double cos_object) const override;
};


class DiskShadowSource: virtual public IrrSource {
protected:
	double relative_semiheight_squared_;
public:
	DiskShadowSource(double relative_semiheight);
	~DiskShadowSource() override = 0;
public:
	double relative_semiheight_squared() const;
	bool shadow(const UnitVec3& direction) const override;
};


class PointAccretorSource: public PointSource, public ConstantAlbedoSource, public DiskShadowSource {
public:
	PointAccretorSource(const Vec3& position, double luminosity, double albedo, double relative_semiheight);
	~PointAccretorSource() override = default;
public:
	PointAccretorSource* clone() const override;
};


class CentralDiskSource: public ElementaryPlainSource, public ConstantAlbedoSource, public DiskShadowSource {
public:
	CentralDiskSource(const Vec3& position, const UnitVec3& plain_normal,
					  double luminosity, double albedo, double relative_semiheight);
	~CentralDiskSource() override = default;
public:
	CentralDiskSource* clone() const override;
};


class IrradiatedStar: public Star {
public:
	using sources_t = std::vector<std::unique_ptr<IrrSource>>;
protected:
	sources_t sources_;
private:
	sources_t cloneSources() const;
public:
	IrradiatedStar(sources_t&& sources, double temp, double radius, unsigned short lod);
	IrradiatedStar(sources_t&& sources, double temp, const RocheLobe& roche_lobe, unsigned short lod);
	IrradiatedStar(const IrradiatedStar& other);
public:
	const sources_t& sources() const;
	void set_sources(sources_t&& value);
	const vald& Qirr() override;
};


#endif //FREDDI_STAR_HPP
