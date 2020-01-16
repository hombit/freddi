#ifndef FREDDI_NS_OPTIONS_HPP
#define FREDDI_NS_OPTIONS_HPP

#include <boost/algorithm/string.hpp> // split

#include <options.hpp>
#include <ns/ns_arguments.hpp>

class NeutronStarOptions: public NeutronStarArguments {
protected:
	static pard fpparamsInitializer(const po::variables_map& vm);
public:
	NeutronStarOptions(const po::variables_map& vm);
	static po::options_description description();
};


class FreddiNeutronStarOptions: public FreddiNeutronStarArguments {
public:
	FreddiNeutronStarOptions(const po::variables_map& vm);
	static po::options_description description();
};

#endif //FREDDI_NS_OPTIONS_HPP
