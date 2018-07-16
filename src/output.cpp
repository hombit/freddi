#include "output.hpp"

#include <sstream>

#include "unit_transfomation.hpp"

constexpr const char FreddiFileOutput::fulldata_header[];

FreddiFileOutput::FreddiFileOutput(const Freddi &freddi_, const boost::program_options::variables_map& vm):
		freddi(&freddi_),
		output(freddi_.args->general->dir + "/" + freddi_.args->general->prefix + ".dat") {
	output << "#t\tMdot\tMdisk\tRhot\tCirrout\tH2R\tTeffout\tTirrout\tQiir2Qvisout\tLx\tmU\tmB\tmV\tmR\tmI\tmJ\t";
	for (int i = 0; i < freddi->args->flux->lambdas.size(); ++i) {
		output << " Fnu" << i;
		for (double j = 0; j < 9 - log10(i + 0.1); ++j) {
			output << "\t";
		}
	}
	output << "\n";
	output << "#days\tg/s\tg\tRsun\tfloat\tfloat\tK\tK\tfloat\terg/s\tmag\tmag\tmag\tmag\tmag\tmag";
	for (int i = 0; i < freddi->args->flux->lambdas.size(); ++i) {
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
		output << "# --rout hadn't been specified, tidal radius " << freddi->args->basic->rout / solar_radius << " Rsun was used"
			   << std::endl;
	}
	output << std::flush;
}

void FreddiFileOutput::dump() {
	auto Nx = freddi->get_state().get_Nx();
	output  << sToDay(freddi->get_state().get_t())
			<< "\t" << freddi->get_state().get_Mdot_in()
			<< "\t" << freddi->get_state().Mdisk()
			<< "\t" << cmToSun(freddi->get_state().get_R()[Nx-1])
			<< "\t" << freddi->get_state().get_Cirr()[Nx-1]
			<< "\t" << freddi->get_state().get_Height()[Nx-1] / freddi->get_state().get_R()[Nx-1]
			<< "\t" << freddi->get_state().get_Tph()[Nx-1]
			<< "\t" << freddi->get_state().get_Tirr()[Nx-1]
			<< "\t" << pow( freddi->get_state().get_Tirr()[Nx-1] / freddi->get_state().get_Tph_vis()[Nx-1], 4. )
			<< "\t" << freddi->get_state().get_Lx()
			<< "\t" << freddi->get_state().mU()
			<< "\t" << freddi->get_state().mB()
			<< "\t" << freddi->get_state().mV()
			<< "\t" << freddi->get_state().mR()
			<< "\t" << freddi->get_state().mI()
			<< "\t" << freddi->get_state().mJ();
	for ( auto &lambda : freddi->args->flux->lambdas ){
		output << "\t" << freddi->get_state().flux(lambda);
	}
	output << std::endl;

	if (freddi->args->general->fulldata) {
		std::ostringstream filename;
		auto i_t = static_cast<int>(std::round(freddi->get_state().get_t() / freddi->args->calc->tau));
		filename << freddi->args->general->dir << "/" << freddi->args->general->prefix << "_" << i_t << ".dat";
		FstreamWithPath output(filename.str());
		output << fulldata_header << sToDay(freddi->get_state().get_t()) << " Mdot_in = " << freddi->get_state().get_Mdot_in() << std::endl;
		for ( int i = 1; i < Nx; ++i ){
			output		<< freddi->get_state().get_h()[i]
				<< "\t" << freddi->get_state().get_R()[i]
				<< "\t" << freddi->get_state().get_F()[i]
				<< "\t" << freddi->get_state().get_Sigma()[i]
				<< "\t" << freddi->get_state().get_Tph()[i]
				<< "\t" << freddi->get_state().get_Tph_vis()[i]
				<< "\t" << freddi->get_state().get_Tirr()[i]
				<< "\t" << freddi->get_state().get_Height()[i]
				<< std::endl;
		}
	}
}
