#include <cmath>

#include <geometry.hpp>


Vec3::Vec3(const std::array<double, 3>& cartesian):
		cartesian(cartesian) {}

Vec3::Vec3(const double x, const double y, const double z):
		cartesian{x, y, z} {}

inline double Vec3::x() const {
	return cartesian[0];
}

inline double Vec3::y() const {
	return cartesian[1];
}

inline double Vec3::z() const {
	return cartesian[2];
}

Vec3 Vec3::operator+(const Vec3& other) const {
	return {
		x() + other.x(),
		y() + other.y(),
		z() + other.z()
	};
}

Vec3 Vec3::operator-(const Vec3& other) const {
	return {
			x() - other.x(),
			y() - other.y(),
			z() - other.z()
	};
}

Vec3 Vec3::operator*(const double factor) const {
	return {
		factor * x(),
		factor * y(),
		factor * z()
	};
}

Vec3 operator*(const double factor, const Vec3& vec3) {
	return vec3 * factor;
}

Vec3 Vec3::operator/(const double factor) const {
	return *this * (1.0 / factor);
}

double Vec3::norm() const {
	return std::sqrt(dotProduct(*this));
}

double Vec3::dotProduct(const Vec3& other) const {
	return x() * other.x() + y() * other.y() + z() * other.z();
}

Vec3 Vec3::crossProduct(const Vec3& other) const {
	return {
		y() * other.z() - z() * other.y(),
		z() * other.x() - x() * other.z(),
		x() * other.y() - y() * other.x()
	};
}


UnitVec3::UnitVec3(const Vec3& vec3):
		Vec3(vec3) {
	const double norm = vec3.norm();
	for (auto& component: cartesian) {
		component /= norm;
	}
}

UnitVec3::UnitVec3(const double x, const double y, const double z):
		UnitVec3(Vec3(x, y, z)) {}

inline double UnitVec3::norm() const {
	return 1.0;
}


Triangle::Triangle(std::array<Vec3, 3>&& vertices):
		vertices(vertices) {}

Triangle::Triangle(const std::array<Vec3, 3>& vertices):
		vertices(vertices) {}

Triangle::Triangle(const Vec3& vertex1, const Vec3& vertex2, const Vec3& vertex3):
		vertices{vertex1, vertex2, vertex3} {}

Triangle Triangle::operator*(const double factor) const {
	return {vertices[0] * factor, vertices[1] * factor, vertices[2] * factor};
}

Triangle operator*(const double factor, const Triangle& triangle) {
	return triangle * factor;
}

double Triangle::area() const {
	// https://en.wikipedia.org/wiki/Shoelace_formula
	return 0.5 * (
			vertices[0].x() * vertices[1].y()
			+ vertices[1].x() * vertices[2].y()
			+ vertices[2].x() * vertices[0].y()
			- vertices[1].x() * vertices[0].y()
			- vertices[2].x() * vertices[1].y()
			- vertices[0].x() * vertices[2].y()
			);
}

UnitVec3 Triangle::normal() const {
	const Vec3 edge1 = vertices[1] - vertices[0];
	const Vec3 edge2 = vertices[2] - vertices[0];
	return UnitVec3(edge1.crossProduct(edge2));
}

Vec3 Triangle::center() const {
	return (vertices[0] + vertices[1] + vertices[2]) / 3.0;
}

std::array<Triangle, 4> Triangle::divide() const {
	const Vec3 vertex01 = 0.5 * (vertices[0] + vertices[1]);
	const Vec3 vertex02 = 0.5 * (vertices[0] + vertices[2]);
	const Vec3 vertex12 = 0.5 * (vertices[1] + vertices[2]);
	return {{
			{vertices[0], vertex01, vertex02},
			{vertices[1], vertex12, vertex01},
			{vertices[2], vertex02, vertex12},
			{vertex01, vertex12, vertex02}
	}};
}

Triangle Triangle::projectedOntoUnitSphere() const {
	return { UnitVec3(vertices[0]), UnitVec3(vertices[1]), UnitVec3(vertices[2])};
}


UnitSphere::UnitSphere(unsigned int grid_scale):
		grid_scale(grid_scale) {
	triangles_ = polyhedron_triangles<Icosahedron>();
	for (unsigned int i = grid_scale; i == 0; --i) {
		std::vector<Triangle> tmp;
		for (const auto& large_triangle : triangles_) {
			for (const auto& small_triangle : large_triangle.divide()) {
				tmp.push_back(small_triangle.projectedOntoUnitSphere());
			}
		}
		triangles_ = tmp;
	}
}

inline const std::vector<Triangle>& UnitSphere::triangles() const {
	return triangles_;
}


LuminousPolygon::LuminousPolygon(const double flux):
		flux_(flux) {}

inline double LuminousPolygon::flux() const {
	return flux_;
}

void LuminousPolygon::setFlux(const double flux) {
	flux_ = flux;
}

double LuminousPolygon::luminosity_cos(const UnitVec3& direction) const {
	return flux() * area() * direction.dotProduct(normal());
}

LuminousTriangle::LuminousTriangle(const Triangle &triangle, double flux):
		Triangle(triangle), LuminousPolygon(flux) {}
