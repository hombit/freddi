#include "output.hpp"

#include <sstream>

#include "unit_transformation.hpp"

constexpr const char FreddiFileOutput::fulldata_header[];

FreddiFileOutput::FreddiFileOutput(FreddiEvolution &freddi_, const boost::program_options::variables_map& vm):
		freddi(&freddi_),
		output(freddi_.state().args.general->dir + "/" + freddi_.state().args.general->prefix + ".dat") {
	output << "#t\tMdot\tMdisk\tRhot\tCirrout\tH2R\tTeffout\tTirrout\tQiir2Qvisout\tLx\tmU\tmB\tmV\tmR\tmI\tmJ\t";
	for (int i = 0; i < freddi_.state().args.flux->lambdas.size(); ++i) {
		output << " Fnu" << i;
		for (double j = 0; j < 9 - log10(i + 0.1); ++j) {
			output << "\t";
		}
	}
	output << "\n";
	output << "#days\tg/s\tg\tRsun\tfloat\tfloat\tK\tK\tfloat\terg/s\tmag\tmag\tmag\tmag\tmag\tmag";
	for (int i = 0; i < freddi_.state().args.flux->lambdas.size(); ++i) {
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
		output << "# --rout hadn't been specified, tidal radius " << freddi_.state().args.basic->rout / solar_radius << " Rsun was used"
			   << std::endl;
	}
	output << std::flush;
}

void FreddiFileOutput::dump() {
	auto state = freddi->state();
	auto Nx = state.Nx();
	output  << sToDay(state.t())
			<< "\t" << state.Mdot_in()
			<< "\t" << state.Mdisk()
			<< "\t" << cmToSun(state.R()[Nx-1])
			<< "\t" << state.Cirr()[Nx-1]
			<< "\t" << state.Height()[Nx-1] / state.R()[Nx-1]
			<< "\t" << state.Tph()[Nx-1]
			<< "\t" << state.Tirr()[Nx-1]
			<< "\t" << pow(state.Tirr()[Nx-1] / state.Tph_vis()[Nx-1], 4. )
			<< "\t" << state.Lx()
			<< "\t" << state.mU()
			<< "\t" << state.mB()
			<< "\t" << state.mV()
			<< "\t" << state.mR()
			<< "\t" << state.mI()
			<< "\t" << state.mJ();
	for ( auto &lambda : state.args.flux->lambdas ){
		output << "\t" << state.flux(lambda);
	}
	output << std::endl;

	if (state.args.general->fulldata) {
		std::ostringstream filename;
		auto i_t = static_cast<int>(std::round(state.t() / state.args.calc->tau));
		filename << state.args.general->dir << "/" << state.args.general->prefix << "_" << i_t << ".dat";
		FstreamWithPath output(filename.str());
		output << fulldata_header << sToDay(state.t()) << " Mdot_in = " << state.Mdot_in() << std::endl;
		for ( int i = 1; i < Nx; ++i ){
			output		<< state.h()[i]
				<< "\t" << state.R()[i]
				<< "\t" << state.F()[i]
				<< "\t" << state.Sigma()[i]
				<< "\t" << state.Tph()[i]
				<< "\t" << state.Tph_vis()[i]
				<< "\t" << state.Tirr()[i]
				<< "\t" << state.Height()[i]
				<< std::endl;
		}
	}
}
