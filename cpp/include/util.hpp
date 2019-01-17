#ifndef FREDDI_UTIL_HPP
#define FREDDI_UTIL_HPP

#include <functional>
#include <vector>


typedef std::vector<double> vecd;

double trapz(const vecd& x, const vecd& y, size_t first, size_t last);
double trapz(const vecd& x, std::function<double (size_t)> f, size_t first, size_t last);

double disk_radial_trapz(const vecd& r, const vecd& y, size_t first, size_t last);
double disk_radial_trapz(const vecd& r, std::function<double (size_t)> f, size_t first, size_t last);


#endif //FREDDI_UTIL_HPP
