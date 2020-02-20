#include <cmath>

#include <geometry.hpp>


Vec3::Vec3(const std::array<double, 3>& cartesian, const double length):
		cartesian_(cartesian), length_(length) {}

Vec3::Vec3(const std::array<double, 3>& cartesian):
		Vec3(cartesian, lengthFromCartesian(cartesian)) {}

Vec3::Vec3(const double x, const double y, const double z):
		Vec3({x, y, z}, lengthFromCartesian({x, y, z})) {}

double Vec3::lengthFromCartesian(const std::array<double, 3>& cartesian) {
	return std::sqrt(cartesian[0] * cartesian[0] + cartesian[1] * cartesian[1] + cartesian[2] * cartesian[2]);
}

double Vec3::x() const {
	return cartesian_[0];
}

double Vec3::y() const {
	return cartesian_[1];
}

double Vec3::z() const {
	return cartesian_[2];
}

double Vec3::length() const {
	return length_;
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

Vec3& Vec3::operator*=(const double factor) {
	for (auto& component : cartesian_) {
		component *= factor;
	}
	length_ *= factor;
	return *this;
}

Vec3 Vec3::operator*(const double factor) const {
	auto vec = *this;
	vec *= factor;
	return vec;
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
		Vec3(vec3 / vec3.length()) {
	length_ = 1.0;
}

UnitVec3::UnitVec3(const double theta, const double phi):
		Vec3(cartesianFromSpherical(theta, phi), 1.0) {}

std::array<double, 3> UnitVec3::cartesianFromSpherical(const double theta, const double phi) {
	const double sin_theta = std::sin(theta);
	return {
			sin_theta * std::cos(phi),
			sin_theta * std::sin(phi),
			std::cos(theta)
	};
}


Triangle::Triangle(std::array<Vec3, 3>&& vertices):
		Triangle(vertices) {}

Triangle::Triangle(const std::array<Vec3, 3>& vertices):
		vertices_(vertices),
		area_(0.5 * (vertices[1] - vertices[0]).crossProduct(vertices[2] - vertices[1]).length()),
		// Use incenter instead?
		// https://en.wikipedia.org/wiki/Incircle_and_excircles_of_a_triangle#Cartesian_coordinates
		center_((vertices[0] + vertices[1] + vertices[2]) / 3.0),
		normal_{(vertices[1] - vertices[0]).crossProduct(vertices[2] - vertices[1])} {}

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

Triangle& Triangle::operator*=(const double factor) {
	for (auto& vertex : vertices_) {
		vertex *= factor;
	}
	area_ *= factor * factor;
	center_ *= factor;
	return *this;
}

Triangle Triangle::operator*(const double factor) const {
	Triangle tr(*this);
	tr *= factor;
	return tr;
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
	return area_;
}

const UnitVec3& Triangle::normal() const {
	return normal_;
}

const Vec3& Triangle::center() const {
	return center_;
}

double Triangle::area_cos(const UnitVec3& uvec3) const {
	const double cos = uvec3.dotProduct(normal());
	if (cos <= 0) {
		return 0;
	}
	return area() * cos;
}

double Triangle::area_cos(const double theta, const double phi) const {
	return area_cos({theta, phi});
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


std::vector<Triangle> unit_sphere_triangles(unsigned short lod) {
	std::vector<Triangle> triangles = polyhedron_triangles<Icosahedron>();
	for (unsigned short i = lod; i != 0; --i) {
		std::vector<Triangle> tmp;
		for (const auto& large_triangle : triangles) {
			for (const auto& small_triangle : large_triangle.divide()) {
				tmp.push_back(small_triangle.projectedOntoUnitSphere());
			}
		}
		triangles = tmp;
	}
	return triangles;
}


std::ostream& operator<<(std::ostream& os, const std::vector<Triangle>& triangles) {
	for (const auto& tr : triangles) {
		 os << tr << "\n";
	}
	return os;
}
