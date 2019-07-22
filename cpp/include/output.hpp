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
};

struct FileOutputLongField {
	std::string name;
	std::string unit;
	std::function<const vecd& ()> func;
};


class BasicFreddiFileOutput {
protected:
	std::shared_ptr<FreddiEvolution> freddi;
	virtual void shortDump();
	virtual void longDump();
	const std::vector<FileOutputShortField> short_fields;
	const std::vector<FileOutputLongField> long_fields;
private:
	FstreamWithPath output;
	std::string fulldata_header;
	std::string initializeFulldataHeader() const;
public:
	BasicFreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi, const boost::program_options::variables_map& vm,
						  const std::vector<FileOutputShortField>&& short_fields, const std::vector<FileOutputLongField>&& long_fields);
	void dump();
	inline std::string path() const { return output.path; }
};

class FreddiFileOutput: public BasicFreddiFileOutput {
private:
	static std::vector<FileOutputShortField> initializeShortFields(const std::shared_ptr<FreddiEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeLongFields(const std::shared_ptr<FreddiEvolution>& freddi);
public:
	FreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi, const boost::program_options::variables_map& vm):
			BasicFreddiFileOutput(freddi, vm, initializeShortFields(freddi), initializeLongFields(freddi)) {}
};


class FreddiNeutronStarFileOutput: public BasicFreddiFileOutput {
private:
	static std::vector<FileOutputShortField> initializeShortFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi);
	static std::vector<FileOutputLongField> initializeLongFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi);
public:
	FreddiNeutronStarFileOutput(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi, const boost::program_options::variables_map& vm):
			BasicFreddiFileOutput(freddi, vm, initializeShortFields(freddi), initializeLongFields(freddi)) {}
};

#endif //FREDDI_FREDDIFILEOUTPUT_H
