#include <fstream> // ifstream

#include <boost/filesystem.hpp> // path

#include "constants.hpp"
#include "passband.hpp"
#include "spectrum.hpp"
#include "unit_transformation.hpp"

std::string Passband::nameFromPath(const std::string& filepath) {
	boost::filesystem::path path = filepath;
	for(; !path.extension().empty(); path = path.stem());
	return path.string();
}

std::vector<Passband::PassbandPoint> Passband::dataFromFile(const std::string& filepath) {
	std::ifstream file;
	file.exceptions(std::ifstream::failbit | std::ifstream::eofbit);
	file.open(filepath);

	std::vector<PassbandPoint> data;
	for (;;) {
		double lambda;
		double transmission;
		try {
			file >> lambda;
			file >> transmission;
		} catch (const std::ifstream::failure& e) {
			break;
		}
		lambda = angstromToCm(lambda);
		data.emplace_back(lambda, transmission);
	}
	return data;
}

vecd Passband::lambdasFromData(const std::vector<PassbandPoint> data) {
	vecd lambdas;
	lambdas.reserve(data.size());
	for (const auto& point : data) {
		lambdas.push_back(point.lambda);
	}
	return lambdas;
}

vecd Passband::transmissionsFromData(const std::vector<PassbandPoint> data) {
	vecd transmissions;
	transmissions.reserve(data.size());
	for (const auto& point : data) {
		transmissions.push_back(point.transmission);
	}
	return transmissions;
}

std::function<double (size_t)> Passband::widthFrequencyIntegrationFunction(const vecd &lambdas, const vecd &transmissions) {
	return [&lambdas, &transmissions](const size_t i) -> double {
		return transmissions[i] * GSL_CONST_CGSM_SPEED_OF_LIGHT / m::pow<2>(lambdas[i]);
	};
}

double Passband::bb_integral(const double temp) const {
	return trapz(lambdas, [this, temp](const size_t i) -> double {
		return transmissions[i] * Spectrum::Planck_lambda(temp, lambdas[i]);
	},
			0,
			data.size() - 1);
}
