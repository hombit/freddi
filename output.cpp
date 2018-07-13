#include "output.hpp"

FreddiFileOutput::FreddiFileOutput(const Freddi &freddi_, const boost::program_options::variables_map& vm):
		freddi(&freddi_),
		output(freddi_.args->general->dir + "/" + freddi_.args->general->prefix + ".dat") {
	output << "#t    Mdot Mdisk Rhot Cirrout H2R   Teffout Tirrout Qiir2Qvisout Lx    mU  mB  mV  mR  mI  mJ ";
	for (int i = 0; i < freddi->args->flux->lambdas.size(); ++i) {
		output << " Fnu" << i;
		for (double j = 0; j < 9 - log10(i + 0.1); ++j) {
			output << " ";
		}
	}
	output << "\n";
	output << "#days g/s  g     Rsun float   float K       K       float        erg/s mag mag mag mag mag mag";
	for (int i = 0; i < freddi->args->flux->lambdas.size(); ++i) {
		output << " erg/s/cm^2/Hz";
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
//	output  << freddi->t / DAY
//			<< "\t" << Mdot_in
//			<< "\t" << Mdisk
//			<< "\t" << R.at(Nx-1) / solar_radius
//			<< "\t" << C_irr
//			<< "\t" << Height.at(Nx-1) / R.at(Nx-1)
//			<< "\t" << Tph.at(Nx-1)
//			<< "\t" << Tirr.at(Nx-1)
//			<< "\t" << pow( Tirr.at(Nx-1) / Tph_vis.at(Nx-1), 4. )
//			<< "\t" << Lx
//			<< "\t" << mU
//			<< "\t" << mB
//			<< "\t" << mV
//			<< "\t" << mR
//			<< "\t" << mI
//			<< "\t" << mJ;
//	for ( auto &lambda : args->flux->lambdas ){
//		output_sum
//			<< "\t" << I_lambda(R, Tph, lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * cosiOverD2;
//	}
//	output_sum      << endl;
}



//	if (freddi->args->general->fulldata){
//		ostringstream filename;
//		filename << args->general->dir << "/" << args->general->prefix << "_" << i_t << ".dat";
//		ofstream output( filename.str() );
//		output << "#h      R  F      Sigma  Teff Tvis Tirr Height" << "\n";
//		output << "#cm^2/s cm dyn*cm g/cm^2 K    K    K    cm" << "\n";
//		output << "# Time = " << t / DAY << " Mdot_in = " << Mdot_in << endl;
//		for ( int i = 1; i < Nx; ++i ){
//			output		<< h.at(i)
//				<< "\t" << R.at(i)
//				<< "\t" << F.at(i)
//				<< "\t" << Sigma.at(i)
//				<< "\t" << Tph.at(i)
//				<< "\t" << Tph_vis.at(i)
//				<< "\t" << Tirr.at(i)
//				<< "\t" << Height.at(i)
//				<< endl;
//		}
//	}
