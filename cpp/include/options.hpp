#ifndef FREDDI_OPTIONS_HPP
#define FREDDI_OPTIONS_HPP


#ifndef INSTALLPATHPREFIX
#define INSTALLPATHPREFIX ""
#endif  // INSTALLPATHPREFIX

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "arguments.hpp"
#include "passband.hpp"
#include "util.hpp"

namespace po = boost::program_options;


class GeneralOptions: public GeneralArguments {
public:
	GeneralOptions(const po::variables_map& vm);
	static po::options_description description();
};


class BasicDiskBinaryOptions: public BasicDiskBinaryArguments {
public:
	BasicDiskBinaryOptions(const po::variables_map& vm);
	static po::options_description description();
	static std::optional<double> rinInitializer(const po::variables_map& vm);
	static std::optional<double> routInitializer(const po::variables_map& vm);
	static std::optional<double> riscoInitializer(const po::variables_map& vm);
};


class DiskStructureOptions: public DiskStructureArguments {
public:
	DiskStructureOptions(const po::variables_map& vm, const BasicDiskBinaryArguments& bdb_args);
	static po::options_description description();
	static pard windparamsInitializer(const po::variables_map& vm);
};


class SelfIrradiationOptions: public SelfIrradiationArguments {
public:
	SelfIrradiationOptions(const po::variables_map& vm, const DiskStructureArguments& dsa_args);
	static po::options_description description();
};


class FluxOptions: public FluxArguments {
protected:
	static vecd lambdasInitializer(const po::variables_map& vm);
	static std::vector<Passband> passbandsInitializer(const po::variables_map& vm);
public:
	FluxOptions(const po::variables_map& vm);
	static po::options_description description();
};


class CalculationOptions: public CalculationArguments {
protected:
	static std::optional<double> tauInitializer(const po::variables_map& vm);
public:
	CalculationOptions(const po::variables_map& vm);
	static po::options_description description();
};


class FreddiOptions: public FreddiArguments {
public:
	FreddiOptions(const po::variables_map& vm);
	static po::options_description description();
};

template <typename T>
std::optional<T> varToOpt(const po::variables_map& vm, const std::string& name) {
	if (vm.count(name) == 0) {
		return {};
	}
	return vm[name].as<T>();
}


template <typename Options>
bool parseOptions(po::variables_map& vm, int ac, char* av[]) {
	const std::string default_config_filename = "freddi.ini";

	const char* xdg_config_home = getenv("XDG_CONFIG_HOME");
	std::string config_home;
	if (xdg_config_home) {
		config_home = xdg_config_home;
	} else {
		config_home = std::string(getenv("HOME")) + "/.config";
	}


	const std::array<std::string, 4> config_file_dirs = {".", config_home, INSTALLPATHPREFIX"/etc", "/etc"};
	std::vector<std::string> config_file_paths;
	for (const auto& dir : config_file_dirs) {
		config_file_paths.push_back(dir + "/" + default_config_filename);
	}

	auto desc = Options::description();
	po::store(po::parse_command_line(ac, av, desc), vm);

	if (vm.count("help") > 0) {
		std::cout << desc << std::endl;
		return false;
	}

	if (vm.count("config") > 0) {
		config_file_paths.insert(config_file_paths.begin(), vm["config"].as<std::string>());
	}

	for (const auto &path : config_file_paths) {
		std::ifstream config(path);
		po::store(po::parse_config_file(config, desc), vm);
	}

	try {
		po::notify(vm);
	} catch(po::error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return false;
	}

	return true;
}

#endif //FREDDI_OPTIONS_HPP
