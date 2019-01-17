#include "util.hpp"

#include <cmath>

double trapz(const vecd& x, const vecd& y, const size_t first, const size_t last) {
	if (first >= last) {
		return 0.;
	}
	double s = y[first] * (x[first + 1] - x[first]) + y[last] * (x[last] - x[last - 1]);
	for (size_t i = first + 1; i <= last - 1; i++) {
		s += y[i] * (x[i + 1] - x[i - 1]);
	}
	return 0.5 * s;
}


double trapz(const vecd& x, std::function<double (size_t)> f, const size_t first, const size_t last) {
	if (first >= last) {
		return 0.;
	}
	double s = f(first) * (x[first + 1] - x[first]) + f(last) * (x[last] - x[last - 1]);
	for (size_t i = first + 1; i <= last - 1; i++) {
		s += f(i) * (x[i + 1] - x[i - 1]);
	}
	return 0.5 * s;
}


double disk_radial_trapz(const vecd& r, const vecd& y, const size_t first, const size_t last) {
	return trapz(r, [&r, &y](const size_t i) -> double { return 2*M_PI * r[i] * y[i]; }, first, last);
}


double disk_radial_trapz(const vecd& r, std::function<double (size_t)> f, const size_t first, const size_t last) {
	return trapz(r, [&r, f](const size_t i) -> double { return 2*M_PI * r[i] * f(i); }, first, last);
}
