#include <cmath>
#include <numeric>

#include <boost/math/special_functions/pow.hpp>

#include <geometry.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testGeometry

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>


namespace m = boost::math;


BOOST_AUTO_TEST_CASE(testVec3_norm) {
	const Vec3 vec1 = {3.0, 4.0, 0.0};
	BOOST_CHECK_CLOSE(vec1.length(), 5.0, 1e-12);

	const Vec3 vec2 = {1.0, 1.0, 1.0};
	BOOST_CHECK_CLOSE(vec2.length(), std::sqrt(3.0), 1e-12);
}

BOOST_AUTO_TEST_CASE(testVec3_dotProduct) {
	const Vec3 vec1 = {1.0, 0.0, 0.0};
	const Vec3 vec2 = {0.0, 1.0, 0.0};
	BOOST_CHECK_EQUAL(vec1.dotProduct(vec2), vec2.dotProduct(vec1));
	BOOST_CHECK_EQUAL(vec1.dotProduct(vec2), 0.0);

	const Vec3 vec3 = {1.0, 2.0, 3.0};
	const Vec3 vec4 = {-3.0, -2.0, -1.0};
	BOOST_CHECK_EQUAL(vec3.dotProduct(vec4), vec4.dotProduct(vec3));
	BOOST_CHECK_EQUAL(vec3.dotProduct(vec4), -10.0);
}

BOOST_AUTO_TEST_CASE(testVec3_norm2_eq_dotProduct) {
	const Vec3 vec = {5.0, 7.0, -3.0};
	BOOST_CHECK_CLOSE(vec.dotProduct(vec), m::pow<2>(vec.length()), 1e-12);
}

BOOST_AUTO_TEST_CASE(testVec3_crossProduct) {
	const Vec3 vec1 = {1.0, 0.0, 0.0};
	const Vec3 vec2 = {0.0, 1.0, 0.0};
	const Vec3 product12 = vec1.crossProduct(vec2);
	const Vec3 minus_product21 = -1.0 * vec2.crossProduct(vec1);
	BOOST_CHECK_EQUAL_COLLECTIONS(product12.cartesian().begin(), product12.cartesian().end(),
								  minus_product21.cartesian().begin(), minus_product21.cartesian().end());
	const auto product12_desired = {0.0, 0.0, 1.0};
	BOOST_CHECK_EQUAL_COLLECTIONS(product12.cartesian().begin(), product12.cartesian().end(),
								  product12_desired.begin(), product12_desired.end());

	const Vec3 vec3 = {1.0, 2.0, 3.0};
	const Vec3 vec4 = {-1.0, -2.0, -3.0};
	const Vec3 product34 = vec3.crossProduct(vec4);
	const Vec3 minus_product43 = -1.0 * vec4.crossProduct(vec3);
	BOOST_CHECK_EQUAL_COLLECTIONS(product34.cartesian().begin(), product34.cartesian().end(),
								  minus_product43.cartesian().begin(), minus_product43.cartesian().end());
	const auto product34_desired = {0.0, 0.0, 0.0};
	BOOST_CHECK_EQUAL_COLLECTIONS(product34.cartesian().begin(), product34.cartesian().end(),
								  product34_desired.begin(), product34_desired.end());
}


BOOST_AUTO_TEST_CASE(testUnitVec3_length) {
	const Vec3 vec1 = {1.0, 2.0, 3.0};
	const UnitVec3 uvec1(vec1);
	BOOST_CHECK_EQUAL(uvec1.length(), 1.0);
	BOOST_CHECK_CLOSE(uvec1.dotProduct(uvec1), 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(testUnitVec3_direction) {
	const Vec3 vec1 = {-1.0, 2.0, -3.0};
	const UnitVec3 uvec1(vec1);
	const Vec3 product = vec1.crossProduct(uvec1);
	const Vec3 product_desired = {0., 0., 0.};
	// May fail due floating point arithmetic
	BOOST_CHECK_EQUAL_COLLECTIONS(product.cartesian().begin(), product.cartesian().end(),
			product_desired.cartesian().begin(), product_desired.cartesian().end());
}

BOOST_AUTO_TEST_CASE(testTriangle_edges) {
	const Triangle tr1(
			{0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0});
	const auto edges1 = tr1.edges();
	BOOST_CHECK_CLOSE(edges1[0].length(), 1.0, 1e-12);
	BOOST_CHECK_CLOSE(edges1[1].length(), std::sqrt(2), 1e-12);
	BOOST_CHECK_CLOSE(edges1[2].length(), 1.0, 1e-12);

	const Triangle tr2(
			{0.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 1.0, std::sqrt(3.0)});
	const auto edges2 = tr2.edges();
	BOOST_CHECK_CLOSE(edges2[0].length(), 2.0, 1e-12);
	BOOST_CHECK_CLOSE(edges2[1].length(), 2.0, 1e-12);
	BOOST_CHECK_CLOSE(edges2[2].length(), 2.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(testTriangle_area) {
	const Triangle tr1(
			{0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0});
	BOOST_CHECK_CLOSE(tr1.area(), 0.5, 1e-12);

	const Triangle tr2(
			{0.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 1.0, std::sqrt(3.0)});
	BOOST_CHECK_CLOSE(tr2.area(), std::sqrt(3.0), 1e-12);
}

BOOST_AUTO_TEST_CASE(testTriangle_normal) {
	const Triangle tr1(
			{0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0});
	const auto normal1 = tr1.normal();
	const auto normal1_desired = {0.0, 0.0, 1.0};
	BOOST_CHECK_EQUAL_COLLECTIONS(normal1.cartesian().begin(), normal1.cartesian().end(),
			normal1_desired.begin(), normal1_desired.end());

	const Triangle tr2(
			{0.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 1.0, std::sqrt(3.0)});
	const auto normal2 = tr2.normal();
	const auto normal2_desired = {1.0, 0.0, 0.0};
	BOOST_CHECK_EQUAL_COLLECTIONS(normal2.cartesian().begin(), normal2.cartesian().end(),
								  normal2_desired.begin(), normal2_desired.end());
}

BOOST_AUTO_TEST_CASE(testIcosahedron_vertices) {
	const auto triangles = polyhedron_triangles<Icosahedron>();
	for (const auto& tr : triangles) {
		for (const auto& vertex : tr.vertices()) {
			BOOST_CHECK_CLOSE(vertex.length(), 1.0, 1e-12);
		}
	}
}

BOOST_AUTO_TEST_CASE(testIcosahedron_edges) {
	const double edge_length = 4.0 / std::sqrt(2.0 * (5.0 + std::sqrt(5.0)));
	const auto triangles = polyhedron_triangles<Icosahedron>();
	for (const auto& tr : triangles) {
		for (const auto& vertex : tr.edges()) {
			BOOST_CHECK_CLOSE(vertex.length(), edge_length, 1e-12);
		}
	}
}

BOOST_AUTO_TEST_CASE(testIcosahedron_triangles) {
	const double edge_length = 4.0 / std::sqrt(2.0 * (5.0 + std::sqrt(5.0)));
	const auto triangles = polyhedron_triangles<Icosahedron>();
	for (const auto& tr : triangles) {
		const auto normal = tr.normal();
		BOOST_CHECK_EQUAL(normal.length(), 1.0);

		const auto center = tr.center();
		BOOST_CHECK_GT(center.length(), 0.0);
		BOOST_CHECK_LT(center.length(), 1.0);

		// Normal should be co-directional with triangles vertices and center
		BOOST_CHECK_GT(normal.dotProduct(center), 0.0);
		for (const auto& vertex : tr.vertices()) {
			BOOST_CHECK_GT(normal.dotProduct(vertex), 0.0);
		}
	}
}


BOOST_AUTO_TEST_CASE(testUnitSphere_0) {
	const auto sph_triangles = unit_sphere_triangles(0);
	const auto icosahedron_triangles = polyhedron_triangles<Icosahedron>();
	BOOST_CHECK_EQUAL_COLLECTIONS(sph_triangles.begin(), sph_triangles.end(),
			icosahedron_triangles.begin(), icosahedron_triangles.end());
}

BOOST_AUTO_TEST_CASE(testUnitSphere_vertices) {
	for (unsigned short grid_scale = 0; grid_scale < 4; ++grid_scale) {
		const auto sph_triangles = unit_sphere_triangles(grid_scale);
		for (const auto& tr : sph_triangles) {
			for (const auto& vertex : tr.vertices()) {
				BOOST_CHECK_CLOSE(vertex.length(), 1.0, 1e-12);
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(testUnitSphere_triangles) {
	for (unsigned short grid_scale = 0; grid_scale < 4; ++grid_scale) {
		const auto sph_triangles = unit_sphere_triangles(grid_scale);
		for (const auto& tr : sph_triangles) {
			const auto normal = tr.normal();
			BOOST_CHECK_EQUAL(normal.length(), 1.0);

			const auto center = tr.center();
			BOOST_CHECK_GT(center.length(), 0.0);
			BOOST_CHECK_LT(center.length(), 1.0);

			// Normal should be co-directional with triangles vertices and center
			BOOST_CHECK_GT(normal.dotProduct(center), 0.0);
			for (const auto& vertex : tr.vertices()) {
				BOOST_CHECK_GT(normal.dotProduct(vertex), 0.0);
			}
		}
	}
}

double triangles_area(const std::vector<Triangle>& triangles) {
	return std::accumulate(triangles.begin(), triangles.end(), 0.0,
			[](double sum, const Triangle& tr) -> double { return sum + tr.area(); });
}

BOOST_AUTO_TEST_CASE(testUnitSphere_area) {
	const double four_pi = 4.0 * M_PI;

	const auto sph0 = unit_sphere_triangles(0);
	double area_smaller = triangles_area(sph0);
	BOOST_CHECK_GT(area_smaller, 0.0);
	BOOST_CHECK_LT(area_smaller, four_pi);

	for (unsigned short grid_scale = 1; grid_scale < 4; ++grid_scale) {
		const auto sph_triangles = unit_sphere_triangles(grid_scale);
		const double area = triangles_area(sph_triangles);
		BOOST_CHECK_GT(area, 0.0);
		BOOST_CHECK_GT(area, area_smaller);
		BOOST_CHECK_LT(area, four_pi);
		area_smaller = area;
	}
}

double triangles_area_cos(const std::vector<Triangle>& triangles, const double theta, const double phi) {
	return std::accumulate(triangles.begin(), triangles.end(), 0.0,
			[theta, phi](double sum, const Triangle& tr) -> double { return sum + tr.area_cos(theta, phi); });
}

BOOST_AUTO_TEST_CASE(testUnitSphere_projected_area) {
	const std::vector<double> thetas = {0.0, M_PI/6.0, 1.0, M_PI/3.0, M_PI/2.0, M_PI};
	const std::vector<double> phis = {0.0, 1.0, M_PI/3.0, 2.0, M_PI, 1.3 * M_PI, 2.0 * M_PI};

	const auto sph0 = unit_sphere_triangles(0);
	std::vector<std::vector<Triangle>> spheres;
	for (unsigned short grid_scale = 1; grid_scale < 4; ++grid_scale) {
		spheres.push_back(unit_sphere_triangles(grid_scale));
	}

	for (const auto& theta : thetas) {
		for (const auto& phi : phis) {
			double area_smaller = triangles_area_cos(sph0, theta, phi);
			BOOST_CHECK_GT(area_smaller, 0.0);
			BOOST_CHECK_LT(area_smaller, M_PI);

			for (const auto& sph : spheres) {
				const double area = triangles_area_cos(sph, theta, phi);
				BOOST_CHECK_GT(area, 0.0);
				BOOST_CHECK_GT(area, area_smaller);
				BOOST_CHECK_LT(area, M_PI);
				area_smaller = area;
			}
		}
	}
}
