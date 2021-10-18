#ifndef FREDDI_NS_OUTPUT_HPP
#define FREDDI_NS_OUTPUT_HPP

#include <output.hpp>
#include <ns/ns_evolution.hpp>

class FreddiNeutronStarFileOutput: public BasicFreddiFileOutput {
public:
	static std::vector<FileOutputShortField> initializeShortFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeDiskStructureFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeStarFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi);
public:
	FreddiNeutronStarFileOutput(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi, const boost::program_options::variables_map& vm):
			BasicFreddiFileOutput(freddi, vm, initializeShortFields(freddi),
					initializeDiskStructureFields(freddi), initializeStarFields(freddi)) {}
};

#endif //FREDDI_NS_OUTPUT_HPP
