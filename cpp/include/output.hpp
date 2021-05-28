#ifndef FREDDI_FREDDIFILEOUTPUT_H
#define FREDDI_FREDDIFILEOUTPUT_H

#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "freddi_evolution.hpp"


// https://stackoverflow.com/a/366969
class FileOrStdoutStream {
private:
	std::unique_ptr<std::ofstream> ofs;
public:
	std::ostream os;
public:
	FileOrStdoutStream();
	FileOrStdoutStream(const std::string& filename);
};


struct FileOutputShortField {
	std::string name;
	std::string unit;
	std::string description;
	std::function<double ()> func;
	FileOutputShortField(const std::string&& name, const std::string&& unit, const std::string&& description, const std::function<double ()>&& func):
			name(name), unit(unit), description(description), func(func) {};
};

struct FileOutputLongField {
	std::string name;
	std::string unit;
	std::string description;
	std::function<const double (size_t)> func;
	FileOutputLongField(const std::string&& name, const std::string&& unit, const std::string&& description, const std::function<const double (size_t)>&& func):
			name(name), unit(unit), description(description), func(func) {};
};


class BasicFreddiFileOutput {
protected:
	std::shared_ptr<FreddiEvolution> freddi;
	const std::vector<FileOutputShortField> short_fields;
	const std::vector<FileOutputLongField> disk_structure_fields;
	const std::vector<FileOutputLongField> star_fields;
protected:
	void shortDump();
	void diskStructureDump();
	void starDump();
private:
	const unsigned short precision;
	FileOrStdoutStream output;
	std::string disk_structure_header;
	std::string star_header;
	static std::string initializeFulldataHeader(const std::vector<FileOutputLongField>& fields);
public:
	BasicFreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi, const boost::program_options::variables_map& vm,
						  std::vector<FileOutputShortField>&& short_fields,
						  std::vector<FileOutputLongField>&& disk_structure_fields,
						  std::vector<FileOutputLongField>&& star_fields);
	void dump();
};

class FreddiFileOutput: public BasicFreddiFileOutput {
public:
	static std::vector<FileOutputShortField> initializeShortFields(const std::shared_ptr<FreddiEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeDiskStructureFields(const std::shared_ptr<FreddiEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeStarFields(const std::shared_ptr<FreddiEvolution>& freddi);
public:
	FreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi, const boost::program_options::variables_map& vm):
			BasicFreddiFileOutput(freddi, vm, initializeShortFields(freddi), initializeDiskStructureFields(freddi),
					initializeStarFields(freddi)) {}
};


#endif //FREDDI_FREDDIFILEOUTPUT_H
