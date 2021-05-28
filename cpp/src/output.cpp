#include "output.hpp"

#include <algorithm>  // max_element
#include <sstream>

#include "unit_transformation.hpp"


FileOrStdoutStream::FileOrStdoutStream():
		os(std::ostream(std::cout.rdbuf())) {}


FileOrStdoutStream::FileOrStdoutStream(const std::string& filename):
		ofs(new std::ofstream(filename)),
		os(ofs->rdbuf()) {}


template <class T>
void outputHeader(std::basic_ostream<char>& os, const T& fields) {
	os << "#" << fields.at(0).name;
	for (size_t i = 1; i < fields.size(); ++i) {
		os << "\t" << fields[i].name;
	}
	os << "\n";

	os << "#" << fields.at(0).unit;
	for (size_t i = 1; i < fields.size(); ++i) {
		os << "\t" << fields[i].unit;
	}
	os << "\n";

	os << "### Columns description\n";
	for (size_t i = 0; i < fields.size(); ++i) {
		const auto& field = fields[i];
		os << "# " << i + 1 << "=" << field.name << " [" << field.unit << "] : " << field.description << "\n";
	}
}


FileOrStdoutStream output_from_args(const FreddiArguments& args) {
	if (args.general->stdout) {
		return FileOrStdoutStream();
	}
	return FileOrStdoutStream(args.general->dir + "/" + args.general->prefix + ".dat");
}


BasicFreddiFileOutput::BasicFreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi,
											 const boost::program_options::variables_map& vm,
											 std::vector<FileOutputShortField>&& short_fields,
											 std::vector<FileOutputLongField>&& disk_structure_fields,
											 std::vector<FileOutputLongField>&& star_fields):
		freddi(freddi),
		precision(freddi->args().general->output_precision),
		output(output_from_args(freddi->args())),
		short_fields(short_fields),
		disk_structure_fields(disk_structure_fields),
		disk_structure_header(initializeFulldataHeader(disk_structure_fields)),
		star_fields(star_fields),
		star_header(initializeFulldataHeader(star_fields)) {
	auto& out = output.os;
	out.precision(precision);

	outputHeader(output.os, short_fields);

	out << "### Parameters\n";
	for (const auto &it : vm) {
		auto &value = it.second.value();
		if (auto v = boost::any_cast<uint32_t>(&value)) {
			out << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<std::string>(&value)) {
			out << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<double>(&value)) {
			out << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<unsigned int>(&value)) {
			out << "# "
				   << it.first.c_str()
				   << "="
				   << *v
				   << "\n";
		} else if (auto v = boost::any_cast<std::vector<double> >(&value)) {
			for (int i = 0; i < v->size(); ++i) {
				out << "# "
					   << it.first.c_str()
					   << "="
					   << v->at(i)
					   << "  # "
					   << i
					   << "\n";
			}
		} else if (auto v = boost::any_cast<std::vector<std::string> >(&value)) {
			for (int i = 0; i < v->size(); ++i) {
				out << "# "
					   << it.first.c_str()
					   << "="
					   << v->at(i)
					   << "  # "
					   << i
					   << "\n";
			}
		} else {
			throw boost::program_options::invalid_option_value(it.first.c_str());
		}
	}

	out << "### Derived values\n"
		<< "# alpha_cold = " << freddi->args().basic->alphacold << "\n"
		<< "# Tidal radius = " << freddi->args().basic->rout / solar_radius << " Rsun\n"
		<< "# ISCO radius = " << freddi->args().basic->risco << " cm\n";

	out << std::flush;
}


std::string BasicFreddiFileOutput::initializeFulldataHeader(const std::vector<FileOutputLongField>& fields) {
	std::ostringstream oss;
	outputHeader(oss, fields);
	return oss.str();
}

void BasicFreddiFileOutput::shortDump() {
	output.os << short_fields[0].func();
	for (size_t i = 1; i < short_fields.size(); ++i) {
		output.os << "\t" << short_fields[i].func();
	}
	output.os << std::endl;
}

void BasicFreddiFileOutput::diskStructureDump() {
	auto filename = (freddi->args().general->dir + "/" + freddi->args().general->prefix
			+ "_" + std::to_string(freddi->i_t()) + ".dat");
	std::ofstream full_output(filename);
	full_output.precision(precision);

	full_output << disk_structure_header
			<< "### t = " << sToDay(freddi->t()) << " days"
			<< std::endl;

	const size_t last = freddi->args().flux->cold_disk ? freddi->Nx() - 1 : freddi->last();

	for ( int i = freddi->first(); i <= last; ++i ){
		full_output << disk_structure_fields.at(0).func(i);
		for (size_t j = 1; j < disk_structure_fields.size(); ++j) {
			full_output << "\t" << disk_structure_fields[j].func(i);
		}
		full_output << "\n";
	}
	full_output << std::flush;
}

void BasicFreddiFileOutput::starDump() {
	auto filename = (freddi->args().general->dir + "/" + freddi->args().general->prefix
					 + "_" + std::to_string(freddi->i_t()) + "_star.dat");
	std::ofstream full_output(filename);
	full_output.precision(precision);

	full_output << star_header
				<< "### t = " << sToDay(freddi->t()) << " days"
				<< std::endl;

	for (size_t i = 0; i < freddi->star().triangles().size(); ++i){
		full_output << star_fields.at(0).func(i);
		for (size_t j = 1; j < star_fields.size(); ++j) {
			full_output << "\t" << star_fields[j].func(i);
		}
		full_output << "\n";
	}
	full_output << std::flush;
}

void BasicFreddiFileOutput::dump() {
	shortDump();

	if (freddi->args().general->fulldata) {
		diskStructureDump();
		if (freddi->args().flux->star) {
			starDump();
		}
	}

}


std::vector<FileOutputShortField> FreddiFileOutput::initializeShortFields(const std::shared_ptr<FreddiEvolution>& freddi) {
	std::vector<FileOutputShortField> fields {
			{"t", "days", "Time moment", [freddi]() {return sToDay(freddi->t());}},
			{"Mdot", "g/s", "Accretion rate onto central object",  [freddi]() {return freddi->Mdot_in();}},
			{"Mdisk", "g", "Mass of the hot disk", [freddi]() {return freddi->Mdisk();}},
			{"Rhot", "Rsun", "Radius of the hot disk", [freddi]() {return cmToSun(freddi->R()[freddi->last()]);}},
			{"Sigmaout", "g/cm^2", "Surface density at the outer radius of the hot disk", [freddi]() {return freddi->Sigma()[freddi->last()];}},
			{"Kirrout", "float", "Irradiation coefficient Kirr at the outer radius of the hot disk", [freddi]() {return freddi->Kirr()[freddi->last()];}},
			{"H2R", "float", "Relative semiheight at the outer radius of the hot disk", [freddi]() {return freddi->Height()[freddi->last()] / freddi->R()[freddi->last()];}},
			{"Teffout", "K", "Effective tempreture at the outer radius of the hot disk", [freddi]() {return freddi->Tph()[freddi->last()];}},
			{"Tirrout", "K", "Irradiation temperature (Qirr / sigma_SB)^1/4 at the outer radius of the hot disk", [freddi]() {return freddi->Tirr()[freddi->last()];}},
			{"Qirr2Qvisout", "float", "Irradiation flux to viscous flux ratio at the outer radius of the hot disk", [freddi]() {return m::pow<4>(freddi->Tirr()[freddi->last()] / freddi->Tph_vis()[freddi->last()]);}},
			{"TphXmax", "keV", "Maximum effective temperature of the disk", [freddi]() {return kToKev(*std::max_element(freddi->Tph_X().begin() + freddi->first(), freddi->Tph_X().begin() + freddi->last() + 1));}},
			{"Lx", "erg/s", "X-ray luminosity of the disk in the given energy range [emin, emax]", [freddi]() {return freddi->Lx();}},
			{"Lbol", "erg/s", "Bolometric luminosity of the disk", [freddi]() {return freddi->Lbol_disk();}},
			{"Fx", "erg/s/cm^2", "X-ray flux of the disk in the given energy range [emin, emax]", [freddi]() {return freddi->Lx() * freddi->angular_dist_disk(freddi->cosi()) / (FOUR_M_PI * m::pow<2>(freddi->distance()));}},
			{"Fbol", "erg/s/cm^2", "Bolometric flux of the disk", [freddi]() {return freddi->Lbol_disk() * freddi->angular_dist_disk(freddi->cosi()) / (FOUR_M_PI * m::pow<2>(freddi->distance()));}},
			{"mU", "mag", "TO BE REMOVED", [freddi]() {return freddi->mU();}},
			{"mB", "mag", "TO BE REMOVED", [freddi]() {return freddi->mB();}},
			{"mV", "mag", "TO BE REMOVED", [freddi]() {return freddi->mV();}},
			{"mR", "mag", "TO BE REMOVED", [freddi]() {return freddi->mR();}},
			{"mI", "mag", "TO BE REMOVED", [freddi]() {return freddi->mI();}},
			{"mJ", "mag", "TO BE REMOVED", [freddi]() {return freddi->mJ();}},
	};
	const bool cold_disk = freddi->args().flux->cold_disk;
	const bool star = freddi->args().flux->star;
	const auto& lambdas = freddi->args().flux->lambdas;
	for (size_t i = 0; i < lambdas.size(); ++i) {
		const double lambda = lambdas[i];
		fields.emplace_back(
				std::string("Fnu") + std::to_string(i),
				"erg/s/cm^2/Hz",
				"Spectral flux density of the hot disk at wavelength of " + std::to_string(cmToAngstrom(lambda)) + " AA",
				[freddi, lambda]() { return freddi->flux(lambda); }
		);
		if (cold_disk) {
			fields.emplace_back(
					std::string("Fnu") + std::to_string(i) + "_cold",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the cold disk at wavelength of " + std::to_string(cmToAngstrom(lambda)) + " AA",
					[freddi, lambda]() { return freddi->flux_region<FreddiState::ColdRegion>(lambda); }
			);
		}
		if (star) {
			fields.emplace_back(
					"Fnu" + std::to_string(i) + "_star",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star at wavelength of " + std::to_string(cmToAngstrom(lambda)) + " AA with respect to an orbital phase",
					[freddi, lambda]() { return freddi->flux_star(lambda); }
			);
			fields.emplace_back(
					"Fnu" + std::to_string(i) + "_star_min",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star at wavelength of " + std::to_string(cmToAngstrom(lambda)) + " AA on the phase of inferior conjunction of the star",
					[freddi, lambda]() { return freddi->flux_star(lambda, 0.0); }
			);
			fields.emplace_back(
					"Fnu" + std::to_string(i) + "_star_max",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star at wavelength of " + std::to_string(cmToAngstrom(lambda)) + " AA on the phase of superior conjunction of the star",
					[freddi, lambda]() { return freddi->flux_star(lambda, M_PI); }
			);
		}
	}
	for (const auto& pb : freddi->args().flux->passbands) {
		fields.emplace_back(
				"Fnu" + pb.name,
				"erg/s/cm^2/Hz",
				"Spectral flux density of the hot disk in passband " + pb.name,
				[freddi, &pb]() { return freddi->flux(pb); }
		);
		if (cold_disk) {
			fields.emplace_back(
					std::string("Fnu") + pb.name + "_cold",
					"Spectral flux density of the cold disk in passband " + pb.name,
					"erg/s/cm^2/Hz",
					[freddi, &pb]() { return freddi->flux_region<FreddiState::ColdRegion>(pb); }
			);
		}
		if (star) {
			fields.emplace_back(
					"Fnu" + pb.name + "_star",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star in passband " + pb.name + " with respect to an orbital phase",
					[freddi, &pb]() { return freddi->flux_star(pb); }
			);
			fields.emplace_back(
					"Fnu" + pb.name + "_star_min",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star in passband " + pb.name + " on the phase of inferior conjunction of the star",
					[freddi, &pb]() { return freddi->flux_star(pb, 0.0); }
			);
			fields.emplace_back(
					"Fnu" + pb.name + "_star_max",
					"erg/s/cm^2/Hz",
					"Spectral flux density of the optical star in passband " + pb.name + " on the phase of superior conjunction of the star",
					[freddi, &pb]() { return freddi->flux_star(pb, M_PI); }
			);
		}
	}
	return fields;
}

std::vector<FileOutputLongField> FreddiFileOutput::initializeDiskStructureFields(const std::shared_ptr<FreddiEvolution>& freddi) {
	return {
			{"h", "cm^2/s", "Keplerian specific angular momentum", [freddi](size_t i) {return freddi->h()[i];}},
			{"R", "cm", "Radius", [freddi](size_t i) {return freddi->R()[i];}},
			{"F", "dyn*cm", "Viscous torque", [freddi](size_t i) {return freddi->F()[i];}},
			{"Sigma", "g/cm^2", "Surface density", [freddi](size_t i) {return freddi->Sigma()[i];}},
			{"Teff", "K", "Effective temperature", [freddi](size_t i) {return freddi->Tph()[i];}},
			{"Tvis", "K", "Viscous temperature (Qvis / sigma_SB)^1/4", [freddi](size_t i) {return freddi->Tph_vis()[i];}},
			{"Tirr", "K", "Irradiation temperature (Qirr / sigma_SB)^1/4", [freddi](size_t i) {return freddi->Tirr()[i];}},
			{"Height", "cm", "Semiheight", [freddi](size_t i) {return freddi->Height()[i];}},
	};
}

std::vector<FileOutputLongField> FreddiFileOutput::initializeStarFields(const std::shared_ptr<FreddiEvolution>& freddi) {
	auto& star = freddi->star();
	const auto& triangles = star.triangles();

	std::vector<FileOutputLongField> fields;

	for (size_t i_vertex = 0; i_vertex < 3; ++i_vertex) {
		const auto name_prefix = "vertex" + std::to_string(i_vertex) + "_";
		const auto desc_suffix = "-coordinate" + std::to_string(i_vertex) + " triangle vertex";
		fields.emplace_back(name_prefix + "x", "cm", "x" + desc_suffix, [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].x();});
		fields.emplace_back(name_prefix + "y", "cm", "y" + desc_suffix, [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].y();});
		fields.emplace_back(name_prefix + "z", "cm", "z" + desc_suffix, [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].z();});
	}

	fields.emplace_back("center_x", "cm", "x-coordinate of triangle center", [triangles](size_t i) {return triangles[i].center().x();});
	fields.emplace_back("center_y", "cm", "y-coordinate of triangle center", [triangles](size_t i) {return triangles[i].center().y();});
	fields.emplace_back("center_z", "cm", "z-coordinate of triangle center", [triangles](size_t i) {return triangles[i].center().z();});

	fields.emplace_back("Tth", "K", "'Thermal' temperature, i.e. temperature of the star without irradiation", [freddi](size_t i) {return freddi->star().Tth()[i];});
	fields.emplace_back("Teff", "K", "'Effective' temperature, i.e. temperature of the star with respect to irradiation", [freddi](size_t i) {return freddi->star().Teff()[i];});

	return fields;
}
