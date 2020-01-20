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
	return integrate([this](size_t i) -> double { return Qirr()[i]; }, direction) / M_PI;
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
		position_(position),
		luminosity_(0.0) {}

double IrrSource::irr_flux(const Vec3& coord, const UnitVec3& normal) const {
	return 0.0;
}

const Vec3& IrrSource::position() const {
	return position_;
}

double IrrSource::luminosity() const {
	return luminosity_;
}

void IrrSource::set_luminosity(double value) {
	luminosity_ = value;
}

double IrrSource::cos_object(const UnitVec3& unit_distance, const UnitVec3& normal) {
	double cos = unit_distance.dotProduct(normal);
	if (cos <= 0.0) {
		return 0.0;
	}
	return cos;
}


double PointSource::irr_flux(const Vec3& coord, const UnitVec3& normal) const {
	const auto vec = coord - position();
	return luminosity() / (FOUR_M_PI * m::pow<2>(vec.r())) * cos_object(UnitVec3(vec), normal);
}


PlainSource::PlainSource(const Vec3& position, const UnitVec3& plain_normal):
		IrrSource(position),
		plain_normal_(plain_normal) {}

double PlainSource::irr_flux(const Vec3& coord, const UnitVec3& normal) const {
	const auto vec = coord - position();
	const auto uvec = UnitVec3(vec);
	const double cos_source = plain_normal().dotProduct(uvec);
	return luminosity() * 2.0 * cos_source / m::pow<2>(vec.r()) * cos_object(uvec, normal);
}

const UnitVec3& PlainSource::plain_normal() const {
	return plain_normal_;
}


IrradiatedStar::IrradiatedStar(std::vector<std::unique_ptr<IrrSource>>&& sources,
		const double temp, const double radius, const unsigned short lod):
		Star(temp, radius, lod),
		sources_(std::move(sources)) {}

const std::vector<std::unique_ptr<IrrSource>>& IrradiatedStar::sources() const {
	return sources_;
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
