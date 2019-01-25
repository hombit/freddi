#include "output.hpp"

#include <sstream>

#include "unit_transformation.hpp"

constexpr const char FreddiFileOutput::fulldata_header[];

FreddiFileOutput::FreddiFileOutput(FreddiEvolution &freddi_, const boost::program_options::variables_map& vm):
		freddi(&freddi_),
		output(freddi_.args().general->dir + "/" + freddi_.args().general->prefix + ".dat") {
	output << "#t\tMdot\tMdisk\tRhot\tCirrout\tH2R\tTeffout\tTirrout\tQiir2Qvisout\tLx\tmU\tmB\tmV\tmR\tmI\tmJ\t";
	for (int i = 0; i < freddi_.args().flux->lambdas.size(); ++i) {
		output << " Fnu" << i;
		for (double j = 0; j < 9 - log10(i + 0.1); ++j) {
			output << "\t";
		}
	}
	output << "\n";
	output << "#days\tg/s\tg\tRsun\tfloat\tfloat\tK\tK\tfloat\terg/s\tmag\tmag\tmag\tmag\tmag\tmag";
	for (int i = 0; i < freddi_.args().flux->lambdas.size(); ++i) {
		output << "\terg/s/cm^2/Hz";
	}
	output << "\n";
	for (const auto &it : vm) {
		auto &value = it.second.value();
		if (auto v = boost::any_cast<uint32_t>(&value)) {
			output << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<std::string>(&value)) {
			output << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<double>(&value)) {
			output << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<unsigned int>(&value)) {
			output << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<std::vector<double> >(&value)) {
			for (int i = 0; i < v->size(); ++i) {
				output << "# "
					   << it.first.c_str()
					   << "="
					   << v->at(i)
					   << "  # "
					   << i
					   << "\n";
			}
		} else {
//			output << "error\n";
			throw boost::program_options::invalid_option_value(it.first.c_str());
		}
	}
	if (vm.count("rout") == 0) {
		output << "# --rout hadn't been specified, tidal radius " << freddi_.args().basic->rout / solar_radius << " Rsun was used"
			   << std::endl;
	}
	output << std::flush;
}

void FreddiFileOutput::dump() {
	auto Nx = freddi->Nx();
	output  << sToDay(freddi->t())
			<< "\t" << freddi->Mdot_in()
			<< "\t" << freddi->Mdisk()
			<< "\t" << cmToSun(freddi->R()[Nx-1])
			<< "\t" << freddi->Cirr()[Nx-1]
			<< "\t" << freddi->Height()[Nx-1] / freddi->R()[Nx-1]
			<< "\t" << freddi->Tph()[Nx-1]
			<< "\t" << freddi->Tirr()[Nx-1]
			<< "\t" << pow(freddi->Tirr()[Nx-1] / freddi->Tph_vis()[Nx-1], 4. )
			<< "\t" << freddi->Lx()
			<< "\t" << freddi->mU()
			<< "\t" << freddi->mB()
			<< "\t" << freddi->mV()
			<< "\t" << freddi->mR()
			<< "\t" << freddi->mI()
			<< "\t" << freddi->mJ();
	for ( auto &lambda : freddi->args().flux->lambdas ){
		output << "\t" << freddi->flux(lambda);
	}
	output << std::endl;

	if (freddi->args().general->fulldata) {
		std::ostringstream filename;
		auto i_t = static_cast<int>(std::round(freddi->t() / freddi->args().calc->tau));
		filename << freddi->args().general->dir << "/" << freddi->args().general->prefix << "_" << i_t << ".dat";
		FstreamWithPath output(filename.str());
		output << fulldata_header << sToDay(freddi->t()) << " Mdot_in = " << freddi->Mdot_in() << std::endl;
		for ( int i = 1; i < Nx; ++i ){
			output		<< freddi->h()[i]
				<< "\t" << freddi->R()[i]
				<< "\t" << freddi->F()[i]
				<< "\t" << freddi->Sigma()[i]
				<< "\t" << freddi->Tph()[i]
				<< "\t" << freddi->Tph_vis()[i]
				<< "\t" << freddi->Tirr()[i]
				<< "\t" << freddi->Height()[i]
				<< std::endl;
		}
	}
}
