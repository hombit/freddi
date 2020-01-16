#ifndef FREDDI_PASSBAND_HPP
#define FREDDI_PASSBAND_HPP

#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "util.hpp"

class Passband {
public:
	struct PassbandPoint {
		double lambda;
		double transmission;
		PassbandPoint(double lambda, double transmission):
				lambda(lambda), transmission(transmission) {}
	};
private:
	static std::string nameFromPath(const std::string& filepath);
	static std::vector<PassbandPoint> dataFromFile(const std::string& filepath);
	static vecd lambdasFromData(const std::vector<PassbandPoint> data);
	static vecd transmissionsFromData(const std::vector<PassbandPoint> data);
	static std::function<double (size_t)> widthFrequencyIntegrationFunction(const vecd& lambdas, const vecd& transmissions);
public:
	const std::string name;
	const std::vector<PassbandPoint> data;
	const vecd lambdas;
	const vecd transmissions;
	const double t_dl;
	const double t_dnu;
public:
	Passband(const std::string& name, const std::vector<PassbandPoint>& data):
			name(name), data(data),
			lambdas(lambdasFromData(data)), transmissions(transmissionsFromData(data)),
			t_dl(trapz(lambdas, transmissions, 0, data.size() - 1)),
			t_dnu(trapz(lambdas, widthFrequencyIntegrationFunction(lambdas, transmissions), 0, data.size() - 1)) {};
	Passband(const std::string& filepath):
			Passband(nameFromPath(filepath), dataFromFile(filepath)) {};
	inline double bb_lambda(double temp) const { return bb_integral(temp) / t_dl; }
	inline double bb_nu(double temp) const { return bb_integral(temp) / t_dnu; }
protected:
	double bb_integral(double temp) const;
};

#endif //FREDDI_PASSBAND_HPP
