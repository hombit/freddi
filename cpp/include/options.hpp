#ifndef FREDDI_OPTIONS_HPP
#define FREDDI_OPTIONS_HPP


#ifndef INSTALLPATHPREFIX
#define INSTALLPATHPREFIX ""
#endif  // INSTALLPATHPREFIX

#include <boost/program_options.hpp>

#include "arguments.hpp"

namespace po = boost::program_options;


class GeneralOptions: public GeneralArguments {
public:
	GeneralOptions(const po::variables_map& vm);
	static po::options_description description();
};


class BasicDiskBinaryOptions: public BasicDiskBinaryArguments {
protected:
	static double rinInitializer(const po::variables_map& vm);
	static double routInitializer(const po::variables_map& vm);
public:
	BasicDiskBinaryOptions(const po::variables_map& vm);
	static po::options_description description();
};


class DiskStructureOptions: public DiskStructureArguments {
protected:
	static double Mdisk0Initializer(const po::variables_map& vm);
	static double Mdot0Initializer(const po::variables_map& vm);
public:
	DiskStructureOptions(const po::variables_map& vm, const BasicDiskBinaryArguments& bdb_args);
	static po::options_description description();
};


class SelfIrradiationOptions: public SelfIrradiationArguments {
public:
	SelfIrradiationOptions(const po::variables_map& vm, const DiskStructureArguments& dsa_args);
	static po::options_description description();
};


class FluxOptions: public FluxArguments {
protected:
	std::vector<double> lambdasInitializer(const po::variables_map& vm) const;
public:
	FluxOptions(const po::variables_map& vm);
	static po::options_description description();
};


class CalculationOptions: public CalculationArguments {
public:
	CalculationOptions(const po::variables_map& vm);
	static po::options_description description();
};


class FreddiOptions: public FreddiArguments {
public:
	FreddiOptions(const po::variables_map& vm);
	static po::options_description description();
};


po::variables_map parseOptions(int ac, char* av[]);


#endif //FREDDI_OPTIONS_HPP
