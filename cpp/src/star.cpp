#include <numeric> // accumulate

#include <constants.hpp>
#include <star.hpp>


Star::Star(const double temp, const double radius, const unsigned short grid_scale):
		triangles_(initializeTriangles(radius, grid_scale)),
		Tth_(triangles_.size(), temp),
		Teff_(),
		luminosity_() {}

std::vector<Triangle> Star::initializeTriangles(const double radius, const unsigned short grid_scale) {
	auto triangles = unit_sphere_triangles(grid_scale);
	for (auto& tr : triangles) {
		tr *= radius;
	}
	return triangles;
}

const std::vector<Triangle>& Star::triangles() const {
	return triangles_;
}

const vecd& Star::Tth() const {
	return Tth_;
}

double Star::integrate(const vecd& values) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area() * values[i];
	}
	return sum;
}

double Star::integrate(std::function<double (size_t)> func) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area() * func(i);
	}
	return sum;
}

double Star::integrate(const vecd& values, const UnitVec3& direction) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area_cos(direction) * values[i];
	}
	return sum;
}

double Star::integrate(std::function<double (size_t)> func, const UnitVec3& direction) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area_cos(direction) * func(i);
	}
	return sum;
}

const vecd& Star::Teff() {
	if (!Teff_) {
		Teff_ = Tth_;
	}
	return *Teff_;
}

double Star::luminosity() {
	if (!luminosity_) {
		luminosity_ = integrate([this](size_t i) -> double {
			return GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT * m::pow<4>(Teff()[i]);
		});
	}
	return *luminosity_;
}
