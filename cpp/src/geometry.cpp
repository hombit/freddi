#include <cmath>

#include <geometry.hpp>


Vec3::Vec3(const std::array<double, 3>& cartesian):
		cartesian_(cartesian) {}

Vec3::Vec3(const double x, const double y, const double z):
		cartesian_{x, y, z} {}

inline double Vec3::x() const {
	return cartesian_[0];
}

inline double Vec3::y() const {
	return cartesian_[1];
}

inline double Vec3::z() const {
	return cartesian_[2];
}

const std::array<double, 3>& Vec3::cartesian() const {
	return cartesian_;
}

bool Vec3::operator==(const Vec3& other) const {
	return (x() == other.x()) && (y() == other.y()) && (z() == other.z());
}

bool Vec3::operator!=(const Vec3& other) const {
	return !(*this == other);
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

std::ostream& operator<<(std::ostream& os, const Vec3& vec3) {
	os << "Vec3["
		<< vec3.x() << ", "
		<< vec3.y() << ", "
		<< vec3.z() << "]";
	return os;
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
		Vec3(vec3.x()/vec3.norm(), vec3.y()/vec3.norm(), vec3.z()/vec3.norm()) {}

UnitVec3::UnitVec3(const double x, const double y, const double z):
		UnitVec3(Vec3(x, y, z)) {}

inline double UnitVec3::norm() const {
	return 1.0;
}


Triangle::Triangle(std::array<Vec3, 3>&& vertices):
		vertices_(vertices) {}

Triangle::Triangle(const std::array<Vec3, 3>& vertices):
		vertices_(vertices) {}

Triangle::Triangle(const Vec3& vertex1, const Vec3& vertex2, const Vec3& vertex3):
		Triangle{{vertex1, vertex2, vertex3}} {}

const std::array<Vec3, 3>& Triangle::vertices() const {
	return vertices_;
}

std::array<Vec3, 3> Triangle::edges() const {
	return {
			vertices_[1] - vertices_[0],
			vertices_[2] - vertices_[1],
			vertices_[0] - vertices_[2]
	};
}

bool Triangle::operator==(const Triangle& other) const {
	return ((vertices_[0] == other.vertices_[0])
		&& (vertices_[1] == other.vertices_[1])
		&& (vertices_[2] == other.vertices_[2]));
}

bool Triangle::operator!=(const Triangle& other) const {
	return !(*this == other);
}

Triangle Triangle::operator*(const double factor) const {
	return {vertices_[0] * factor, vertices_[1] * factor, vertices_[2] * factor};
}

Triangle operator*(const double factor, const Triangle& triangle) {
	return triangle * factor;
}

std::ostream& operator<<(std::ostream& os, const Triangle& triangle) {
	os << "Triangle["
	   << triangle.vertices()[0] << ", "
	   << triangle.vertices()[1] << ", "
	   << triangle.vertices()[2] << "]";
	return os;
}

double Triangle::area() const {
	const Vec3 edge1 = vertices_[1] - vertices_[0];
	const Vec3 edge2 = vertices_[2] - vertices_[0];
	return 0.5 * edge1.crossProduct(edge2).norm();
}

UnitVec3 Triangle::normal() const {
	const Vec3 edge1 = vertices_[1] - vertices_[0];
	const Vec3 edge2 = vertices_[2] - vertices_[0];
	return UnitVec3(edge1.crossProduct(edge2));
}

// Use incenter instead?
// https://en.wikipedia.org/wiki/Incircle_and_excircles_of_a_triangle#Cartesian_coordinates
Vec3 Triangle::center() const {
	return (vertices_[0] + vertices_[1] + vertices_[2]) / 3.0;
}

std::array<Triangle, 4> Triangle::divide() const {
	const Vec3 vertex01 = 0.5 * (vertices_[0] + vertices_[1]);
	const Vec3 vertex02 = 0.5 * (vertices_[0] + vertices_[2]);
	const Vec3 vertex12 = 0.5 * (vertices_[1] + vertices_[2]);
	return {{
			{vertices_[0], vertex01, vertex02},
			{vertices_[1], vertex12, vertex01},
			{vertices_[2], vertex02, vertex12},
			{vertex01, vertex12, vertex02}
	}};
}

Triangle Triangle::projectedOntoUnitSphere() const {
	return {UnitVec3(vertices_[0]), UnitVec3(vertices_[1]), UnitVec3(vertices_[2])};
}


UnitSphere::UnitSphere(unsigned short grid_scale):
		grid_scale(grid_scale) {
	triangles_ = polyhedron_triangles<Icosahedron>();
	for (unsigned short i = grid_scale; i != 0; --i) {
		std::vector<Triangle> tmp;
		for (const auto& large_triangle : triangles_) {
			for (const auto& small_triangle : large_triangle.divide()) {
				tmp.push_back(small_triangle.projectedOntoUnitSphere());
			}
		}
		triangles_ = tmp;
	}
}

const std::vector<Triangle>& UnitSphere::triangles() const {
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
	const double cos = direction.dotProduct(normal());
	if (cos <= 0) {
		return 0;
	}
	return flux() * area() * cos;
}

LuminousTriangle::LuminousTriangle(const Triangle &triangle, double flux):
		Triangle(triangle), LuminousPolygon(flux) {}
