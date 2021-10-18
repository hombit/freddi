#ifndef FREDDI_GEOMETRY_HPP
#define FREDDI_GEOMETRY_HPP

#include <array>
#include <ostream>
#include <vector>

#include <arguments.hpp>

class Vec3 {
protected:
	std::array<double, 3> cartesian_;
	double length_;
private:
	static double lengthFromCartesian(const std::array<double, 3>& cartesian);
protected:
	Vec3(const std::array<double, 3>& cartesian, double length);
public:
	explicit Vec3(const std::array<double, 3>& cartesian);
	Vec3(double x, double y, double z);
	Vec3(const Vec3& other) = default;
public:
	double x() const;
	double y() const;
	double z() const;
	double length() const;
	const std::array<double, 3>& cartesian() const;
	bool operator==(const Vec3& other) const;
	bool operator!=(const Vec3& other) const;
	Vec3 operator+(const Vec3& other) const;
	Vec3 operator-(const Vec3& other) const;
	Vec3& operator*=(double factor);
	Vec3 operator*(double factor) const;
	Vec3 operator/(double factor) const;
	double dotProduct(const Vec3& other) const;
	Vec3 crossProduct(const Vec3& other) const;
};

Vec3 operator*(double factor, const Vec3& vec3);
std::ostream& operator<<(std::ostream& os, const Vec3& vec3);


class UnitVec3: public Vec3 {
private:
	static std::array<double, 3> cartesianFromSpherical(double theta, double phi);
public:
	explicit UnitVec3(const Vec3& vec3);
	UnitVec3(double theta, double phi);
};


class Triangle {
protected:
	std::array<Vec3, 3> vertices_;
	double area_;
	Vec3 center_;
	UnitVec3 normal_;
public:
	Triangle(std::array<Vec3, 3>&& vertices);
	Triangle(const std::array<Vec3, 3>& vertices);
	Triangle(const Vec3& vertex1, const Vec3& vertex2, const Vec3& vertex3);
public:
	const std::array<Vec3, 3>& vertices() const;
	std::array<Vec3, 3> edges() const;
	bool operator==(const Triangle& other) const;
	bool operator!=(const Triangle& other) const;
	Triangle& operator*=(double fector);
	Triangle operator*(double factor) const;
	double area() const;
	const UnitVec3& normal() const;
	const Vec3& center() const;
	double area_cos(const UnitVec3& uvec3) const;
	double area_cos(double theta, double phi) const;
	std::array<Triangle, 4> divide() const;
	Triangle projectedOntoUnitSphere() const;
};

Triangle operator*(double factor, const Triangle& triangle);
std::ostream& operator<<(std::ostream& os, const Triangle& triangle);


template <size_t VertexCount, size_t TriangleCount>
class TriangularPolyhedron {
public:
	static const constexpr size_t vertex_count = VertexCount;
	static const constexpr size_t triangle_count = TriangleCount;
	const std::array<std::array<double, 3>, VertexCount> vertices;
	const std::array<std::array<size_t, 3>, TriangleCount> triangle_indices;
};

template <class Polyhedron>
std::vector<Triangle> polyhedron_triangles() {
	std::vector<Triangle> res;
	res.reserve(Polyhedron::triangle_count);
	for (const auto& indices : Polyhedron::triangle_indices) {
		res.emplace_back(
				Vec3(Polyhedron::vertices[indices[0]]),
				Vec3(Polyhedron::vertices[indices[1]]),
				Vec3(Polyhedron::vertices[indices[2]]));
	}
	return res;
}


// http://www.opengl.org.ru/docs/pg/0208.html
class Icosahedron: public TriangularPolyhedron<12, 20> {
public:
	static const constexpr double x = .525731112119133606;
	static const constexpr double z = .850650808352039932;
	static const constexpr std::array<std::array<double, 3>, vertex_count> vertices = {{
		{-x, 0.0, z}, {x, 0.0, z}, {-x, 0.0, -z}, {x, 0.0, -z},
		{0.0, z, x}, {0.0, z, -x}, {0.0, -z, x}, {0.0, -z, -x},
		{z, x, 0.0}, {-z, x, 0.0}, {z, -x, 0.0}, {-z, -x, 0.0}
	}};
	static const constexpr std::array<std::array<size_t, 3>, triangle_count> triangle_indices = {{
		{0,1,4}, {0,4,9}, {9,4,5}, {4,8,5}, {4,1,8},
		{8,1,10}, {8,10,3}, {5,8,3}, {5,3,2}, {2,3,7},
		{7,3,10}, {7,10,6}, {7,6,11}, {11,6,0}, {0,6,1},
		{6,10,1}, {9,11,0}, {9,2,11}, {9,5,2}, {7,11,2}
	}};
};


std::vector<Triangle> unit_sphere_triangles(unsigned short lod);


std::ostream& operator<<(std::ostream& os, const std::vector<Triangle>& triangles);

#endif //FREDDI_GEOMETRY_HPP
