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


class FreddiFileOutput {
protected:
	struct ShortField {
		std::string name;
		std::string unit;
		std::function<double ()> func;
	};
	struct LongField {
		std::string name;
		std::string unit;
		std::function<const vecd& ()> func;
	};

	FreddiEvolution* freddi;
	virtual void shortDump();
	virtual void longDump();
	const std::vector<ShortField> short_fields;
	const std::vector<LongField> long_fields;
	virtual std::vector<ShortField> initializeShortFields() const;
	virtual std::vector<LongField> initializeLongFields() const;

private:
	FstreamWithPath output;
	std::string fulldata_header;
	std::string initializeFulldataHeader() const;
public:
	FreddiFileOutput(FreddiEvolution& freddi, const boost::program_options::variables_map& vm);
	void dump();
	inline std::string path() const { return output.path; }
};


class FreddiNeutronStarFileOutput: public FreddiFileOutput {

};


#endif //FREDDI_FREDDIFILEOUTPUT_H
