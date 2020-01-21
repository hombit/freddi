#include <numeric> // accumulate

#include <constants.hpp>
#include <spectrum.hpp>
#include <star.hpp>
#include <util.hpp>


Star::Star(const double temp, const double radius, const unsigned short lod):
		triangles_(initializeTriangles(radius, lod)),
		Tth_(temp, triangles_.size()),
		irr_() {}

std::vector<Triangle> Star::initializeTriangles(const double radius, const unsigned short grid_scale) {
	auto triangles = unit_sphere_triangles(grid_scale);
	for (auto& tr : triangles) {
		tr *= radius;
	}
	return triangles;
}

void Star::invalidate_irradiated_properties() {
	irr_ = IrradiatedProperties();
}

const std::vector<Triangle>& Star::triangles() const {
	return triangles_;
}

const vald& Star::Tth() const {
	return Tth_;
}

double Star::integrate(std::function<double (size_t)>&& func) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area() * func(i);
	}
	return sum;
}

double Star::integrate(std::function<double (size_t)>&& func, const UnitVec3& direction) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area_cos(direction) * func(i);
	}
	return sum;
}

const vald& Star::Qirr() {
	if (!irr_.Qirr) {
		irr_.Qirr = vald(0.0, triangles().size());
	}
	return *irr_.Qirr;
}

const vald& Star::Teff() {
	if (!irr_.Teff) {
		irr_.Teff = std::pow(Qirr() / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT + m::pow<4>(Tth()), 0.25);
	}
	return *irr_.Teff;
}

double Star::luminosity() {
	if (!irr_.luminosity) {
		irr_.luminosity = integrate([this](size_t i) -> double {
			return GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * m::pow<4>(Teff()[i]);
		});
	}
	return *irr_.luminosity;
}

double Star::dLdOmega(const UnitVec3& direction) {
	const double integral = integrate([this](size_t i) -> double {
			return m::pow<4>(Teff()[i]);
		}, direction);
	return GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT / M_PI * integral;
}

double Star::dLdOmega(const UnitVec3& direction, const double lambda) {
	return integrate([this, lambda](size_t i) -> double {
			return Spectrum::Planck_lambda(Teff()[i], lambda) * m::pow<2>(lambda) / GSL_CONST_CGSM_SPEED_OF_LIGHT;
		},
			direction);
}

double Star::dLdOmega(const UnitVec3& direction, const Passband& passband) {
	return integrate([this, &passband](size_t i) -> double { return passband.bb_nu(Teff()[i]); },
			direction);
}

IrrSource::IrrSource(const Vec3& position):
		position_(position) {}

IrrSource::~IrrSource()	{}

double IrrSource::irr_flux(const Vec3& coord, const UnitVec3& normal) const {
	const auto distance = coord - position();
	const UnitVec3 direction(distance);
	if (shadow(direction)) {
		return 0;
	}
	return irr_dLdOmega(direction) / m::pow<2>(distance.r()) * cos_object(direction, normal);
}

const Vec3& IrrSource::position() const {
	return position_;
}

double IrrSource::cos_object(const UnitVec3& unit_distance, const UnitVec3& normal) {
	double cos = unit_distance.dotProduct(normal);
	if (cos <= 0.0) {
		return 0.0;
	}
	return cos;
}


PointLikeSource::PointLikeSource(const Vec3& position, const double luminosity):
		IrrSource(position),
		luminosity_(luminosity) {}

PointLikeSource::~PointLikeSource() {}

double PointLikeSource::luminosity() const {
	return luminosity_;
}


double PointSource::irr_dLdOmega(const UnitVec3& direction) const {
	return luminosity() / FOUR_M_PI;
}

PointSource::~PointSource() {}


ElementaryPlainSource::ElementaryPlainSource(const Vec3& position, const UnitVec3& plain_normal, const double luminosity):
		PointLikeSource(position, luminosity),
		plain_normal_(plain_normal) {}

ElementaryPlainSource::~ElementaryPlainSource() {}

const UnitVec3& ElementaryPlainSource::plain_normal() const {
	return plain_normal_;
}

double ElementaryPlainSource::irr_dLdOmega(const UnitVec3& direction) const {
	const double cos_source = std::abs(plain_normal().dotProduct(direction));
	return luminosity() * 2.0 * cos_source;
}


DiskShadowSource::DiskShadowSource(const Vec3& position, double relative_semiheight):
		IrrSource(position),
		relative_semiheight_squared_(m::pow<2>(relative_semiheight)) {}

DiskShadowSource::~DiskShadowSource() {}

bool DiskShadowSource::shadow(const UnitVec3& direction) const {
	const double rho2 = m::pow<2>(direction.x()) + m::pow<2>(direction.y());
	return m::pow<2>(direction.z()) / rho2 < relative_semiheight_squared();
}

double DiskShadowSource::relative_semiheight_squared() const {
	return relative_semiheight_squared_;
}


PointAccretorSource::PointAccretorSource(const Vec3& position, double luminosity, double relative_semiheight):
		IrrSource(position),
		PointSource(position, luminosity),
		DiskShadowSource(position, relative_semiheight) {}

PointAccretorSource* PointAccretorSource::clone() const {
	return new PointAccretorSource(*this);
}

CentralDiskSource::CentralDiskSource(const Vec3& position, const UnitVec3& plain_normal,
									 const double luminosity, const double relative_semiheight):
		IrrSource(position),
		ElementaryPlainSource(position, plain_normal, luminosity),
		DiskShadowSource(position, relative_semiheight) {}

CentralDiskSource* CentralDiskSource::clone() const {
	return new CentralDiskSource(*this);
}

IrradiatedStar::IrradiatedStar(sources_t&& sources, const double temp, const double radius, const unsigned short lod):
		Star(temp, radius, lod),
		sources_(std::move(sources)) {}

IrradiatedStar::IrradiatedStar(const IrradiatedStar& other):
		Star(other),
		sources_(cloneSources()) {}

IrradiatedStar::sources_t IrradiatedStar::cloneSources() const {
	sources_t result;
	for (const auto& source : sources()) {
		result.emplace_back(source->clone());
	}
	return result;
}

const IrradiatedStar::sources_t& IrradiatedStar::sources() const {
	return sources_;
}

void IrradiatedStar::set_sources(IrradiatedStar::sources_t&& value) {
	invalidate_irradiated_properties();
	sources_ = std::move(value);
}

const vald& IrradiatedStar::Qirr() {
	if (!irr_.Qirr) {
		vald x(0.0, triangles().size());
		for (size_t i = 0; i < triangles().size(); ++i) {
			for (const auto& source : sources()) {
				x[i] += source->irr_flux(triangles_[i].center(), triangles_[i].normal());
			}
		}
		irr_.Qirr = std::move(x);
	}
	return *irr_.Qirr;
}
