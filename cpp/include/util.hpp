#ifndef FREDDI_UTIL_HPP
#define FREDDI_UTIL_HPP

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#include <boost/math/special_functions/pow.hpp>


namespace m = boost::math;


typedef std::valarray<double> vald;
typedef std::vector<double> vecd;
typedef std::map<std::string, double> pard;

double trapz(const vecd& x, const vecd& y, size_t first, size_t last);
double trapz(const vecd& x, const std::function<double (size_t)>& f, size_t first, size_t last);

double disk_radial_trapz(const vecd& r, const vecd& y, size_t first, size_t last);
double disk_radial_trapz(const vecd& r, const std::function<double (size_t)>& f, size_t first, size_t last);

double simps(const vecd& x, const vecd& y, size_t first, size_t last);
double simps(const vecd& x, const std::function<double (size_t)>& f, size_t first, size_t last);

#endif //FREDDI_UTIL_HPP
