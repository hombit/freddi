#ifndef FREDDI_FREDDIFILEOUTPUT_H
#define FREDDI_FREDDIFILEOUTPUT_H

#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "freddi_evolution.hpp"


class FstreamWithPath: public std::ofstream {
public:
	const std::string path;
public:
	FstreamWithPath(const std::string& path):
			std::ofstream(path),
			path(path) {}
};

struct FileOutputShortField {
	std::string name;
	std::string unit;
	std::function<double ()> func;
	FileOutputShortField(const std::string&& name, const std::string&& unit, const std::function<double ()>&& func):
			name(name), unit(unit), func(func) {};
};

struct FileOutputLongField {
	std::string name;
	std::string unit;
	std::function<const double (size_t)> func;
	FileOutputLongField(const std::string&& name, const std::string&& unit, const std::function<const double (size_t)>&& func):
			name(name), unit(unit), func(func) {};
};


class BasicFreddiFileOutput {
protected:
	constexpr static const int precision = 12;
	std::shared_ptr<FreddiEvolution> freddi;
	void shortDump();
	void diskStructureDump();
	void starDump();
	const std::vector<FileOutputShortField> short_fields;
	const std::vector<FileOutputLongField> disk_structure_fields;
	const std::vector<FileOutputLongField> star_fields;
private:
	FstreamWithPath output;
	std::string disk_structure_header;
	std::string star_header;
	static std::string initializeFulldataHeader(const std::vector<FileOutputLongField>& fields);
public:
	BasicFreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi, const boost::program_options::variables_map& vm,
						  std::vector<FileOutputShortField>&& short_fields,
						  std::vector<FileOutputLongField>&& disk_structure_fields,
						  std::vector<FileOutputLongField>&& star_fields);
	void dump();
	inline std::string path() const { return output.path; }
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
