#ifndef FREDDI_FREDDIFILEOUTPUT_H
#define FREDDI_FREDDIFILEOUTPUT_H

#include <fstream>
#include <string>

#include <boost/program_options.hpp>

#include "freddi.hpp"


class FstreamWithPath: public std::ofstream {
public:
	const std::string path;
public:
	FstreamWithPath(const std::string& path):
			std::ofstream(path),
			path(path) {}
};


class FreddiFileOutput {
private:
	const Freddi* freddi;
	FstreamWithPath output;
public:
	FreddiFileOutput(const Freddi& freddi, const boost::program_options::variables_map& vm);
	void dump();
	inline std::string path() const { return output.path; }
};


#endif //FREDDI_FREDDIFILEOUTPUT_H
