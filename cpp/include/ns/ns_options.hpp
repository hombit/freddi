#ifndef FREDDI_NS_OPTIONS_HPP
#define FREDDI_NS_OPTIONS_HPP

#include <options.hpp>
#include <ns/ns_arguments.hpp>

class NeutronStarOptions: public NeutronStarArguments {
protected:
	static pard fpparamsInitializer(const po::variables_map& vm);
	static pard kappatparamsInitalizer(const po::variables_map& vm);
public:
	NeutronStarOptions(const po::variables_map& vm);
	static po::options_description description();
};

class NeutronStarBasicDiskBinaryOptions: public NeutronStarBasicDiskBinaryArguments {
public:
	NeutronStarBasicDiskBinaryOptions(const po::variables_map& vm, const NeutronStarArguments& ns_args);
	static po::options_description description();
};

class NeutronStarSelfIrradiationOptions: public NeutronStarSelfIrradiationArguments {
public:
	NeutronStarSelfIrradiationOptions(const po::variables_map& vm, const DiskStructureArguments& dsa_args);
	static po::options_description description();
};

class FreddiNeutronStarOptions: public FreddiNeutronStarArguments {
public:
	FreddiNeutronStarOptions(const po::variables_map& vm);
	static po::options_description description();
};

#endif //FREDDI_NS_OPTIONS_HPP
