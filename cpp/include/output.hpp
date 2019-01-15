#ifndef FREDDI_FREDDIFILEOUTPUT_H
#define FREDDI_FREDDIFILEOUTPUT_H

#include <fstream>
#include <string>

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
	constexpr static const char fulldata_header[] = "#h\tR\tF\tSigma\tTeff\tTvis\tTirr\tHeight\n#cm^2/s\tcm\tdyn*cm\tg/cm^2\tK\tK\tK\tcm\n# Time = ";
private:
	FreddiEvolution* freddi;
	FstreamWithPath output;
public:
	FreddiFileOutput(FreddiEvolution& freddi, const boost::program_options::variables_map& vm);
	void dump();
	inline std::string path() const { return output.path; }
};


#endif //FREDDI_FREDDIFILEOUTPUT_H
