#include <numeric> // accumulate

#include <constants.hpp>
#include <star.hpp>


Star::Star(const double temp, const double radius, const unsigned short grid_scale):
		triangles_(initializeTriangles(radius, grid_scale)),
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

double Star::integrate(std::function<double (size_t)> func) const {
	double sum = 0;
	for (size_t i = 0; i < triangles().size(); ++i) {
		sum += triangles()[i].area() * func(i);
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

const vald& Star::Tirr() {
	if (!irr_.Tirr) {
		irr_.Tirr = vald(0.0, triangles().size());
	}
	return *irr_.Tirr;
}

const vald& Star::Teff() {
	if (!irr_.Teff) {
		irr_.Teff = std::pow(m::pow<4>(Tth()) + m::pow<4>(Tirr()), 0.25);
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


FreddiStar::FreddiStar(const FreddiArguments& args):
		Star(args.basic->Topt, args.basic->Ropt, 3),
		args(args) {}
