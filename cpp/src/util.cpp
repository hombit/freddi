#include "util.hpp"

#include <cmath>

#include <boost/math/special_functions/pow.hpp>

namespace bm = boost::math;

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


double trapz(const vecd& x, const std::function<double (size_t)>& f, const size_t first, const size_t last) {
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


double disk_radial_trapz(const vecd& r, const std::function<double (size_t)>& f, const size_t first, const size_t last) {
	return trapz(r, [&r, f](const size_t i) -> double { return 2*M_PI * r[i] * f(i); }, first, last);
}


double simps(const vecd& x, const vecd& y, const size_t first, const size_t last) {
	const size_t N = last - first + 1;
	switch (N) {
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return 0.5 * (y[first] + y[last]) * (x[last] - x[first]);
	}
	if (N % static_cast<size_t>(2) == static_cast<size_t>(0)) {
		return simps(x, y, first + 1, last) + simps(x, y, first, first + 1);
	}

	double delta = (x[first+2] - x[first]);
	double s = y[first] * delta * (2. - (x[first+2] - x[first+1]) / (x[first+1] - x[first]))
			+ y[first+1] * bm::pow<3>(delta) / ((x[first+2] - x[first+1]) * (x[first+1] - x[first]))
			+ y[last] * (2. - (x[last-1] - x[last]) / (x[last] - x[last] - 1));
	for (size_t i = first + 2; i <= last - 2; i += 2) {
		delta = x[i+2] - x[i];
		s += y[i] * (delta * (2. - (x[i+2] - x[i+1]) / (x[i+1] - x[i]))
				+ (x[i] - x[i-2]) * (2. - (x[i-1] - x[i-2]) / (x[i] - x[i-1])))
			+ y[i+1] * bm::pow<3>(delta) / ((x[i+2] - x[i+1]) * (x[i+1] - x[i]));
	}
	s /= 6.;
	return s;
}


double simps(const vecd& x, const std::function<double (size_t)>& f, const size_t first, const size_t last) {
	const size_t N = last - first + 1;
	switch (N) {
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return 0.5 * (f(first) + f(last)) * (x[last] - x[first]);
	}
	if (N % static_cast<size_t>(2) == static_cast<size_t>(0)) {
		return simps(x, f, first + 1, last) + simps(x, f, first, first + 1);
	}

	double delta = (x[first+2] - x[first]);
	double s = f(first) * delta * (2. - (x[first+2] - x[first+1]) / (x[first+1] - x[first]))
			   + f(first+1) * bm::pow<3>(delta) / ((x[first+2] - x[first+1]) * (x[first+1] - x[first]))
			   + f(last) * (2. - (x[last-1] - x[last]) / (x[last] - x[last] - 1));
	for (size_t i = first + 2; i <= last - 2; i += 2) {
		delta = x[i+2] - x[i];
		s += f(i) * (delta * (2. - (x[i+2] - x[i+1]) / (x[i+1] - x[i]))
					 + (x[i] - x[i-2]) * (2. - (x[i-1] - x[i-2]) / (x[i] - x[i-1])))
			 + f(i+1) * bm::pow<3>(delta) / ((x[i+2] - x[i+1]) * (x[i+1] - x[i]));
	}
	s /= 6.;
	return s;
}
