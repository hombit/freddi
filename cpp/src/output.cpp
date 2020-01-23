#include "output.hpp"

#include <sstream>

#include "unit_transformation.hpp"

BasicFreddiFileOutput::BasicFreddiFileOutput(const std::shared_ptr<FreddiEvolution>& freddi,
											 const boost::program_options::variables_map& vm,
											 std::vector<FileOutputShortField>&& short_fields,
											 std::vector<FileOutputLongField>&& disk_structure_fields,
											 std::vector<FileOutputLongField>&& star_fields):
		freddi(freddi),
		output(freddi->args().general->dir + "/" + freddi->args().general->prefix + ".dat"),
		short_fields(short_fields),
		disk_structure_fields(disk_structure_fields),
		disk_structure_header(initializeFulldataHeader(disk_structure_fields)),
		star_fields(star_fields),
		star_header(initializeFulldataHeader(star_fields)) {
	output << "#" << short_fields.at(0).name;
	for (size_t i = 1; i < short_fields.size(); ++i) {
		output << "\t" << short_fields[i].name;
	}
	output << "\n";

	output << "#" << short_fields.at(0).unit;
	for (size_t i = 1; i < short_fields.size(); ++i) {
		output << "\t" << short_fields[i].unit;
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
		} else if (auto v = boost::any_cast<std::vector<std::string> >(&value)) {
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
		output << "# --rout hadn't been specified, tidal radius " << freddi->args().basic->rout / solar_radius << " Rsun was used"
			   << std::endl;
	}
	output << std::flush;
}


std::string BasicFreddiFileOutput::initializeFulldataHeader(const std::vector<FileOutputLongField>& fields) {
	std::string s;

	s += "#" + fields.at(0).name;
	for (size_t i = 1; i < fields.size(); ++i) {
		s += "\t" + fields[i].name;
	}
	s += "\n";

	s += "#" + fields.at(0).unit;
	for (size_t i = 1; i < fields.size(); ++i) {
		 s += "\t" + fields[i].unit;
	}
	s += "\n";

	return s;
}

void BasicFreddiFileOutput::shortDump() {
	output << short_fields[0].func();
	for (size_t i = 1; i < short_fields.size(); ++i) {
		output << "\t" << short_fields[i].func();
	}
	output << std::endl;
}

void BasicFreddiFileOutput::diskStructureDump() {
	auto filename = (freddi->args().general->dir + "/" + freddi->args().general->prefix
			+ "_" + std::to_string(freddi->i_t()) + ".dat");
	FstreamWithPath full_output(filename);

	full_output << disk_structure_header
			<< "# Time = " << sToDay(freddi->t())
			<< " Mdot_in = " << freddi->Mdot_in()
			<< std::endl;

	for ( int i = freddi->first(); i <= freddi->last(); ++i ){
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
	FstreamWithPath full_output(filename);

	full_output << star_header
				<< "# Time = " << sToDay(freddi->t())
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
			{"t", "days", [freddi]() {return sToDay(freddi->t());}},
			{"Mdot", "g/s",  [freddi]() {return freddi->Mdot_in();}},
			{"Mdisk", "g", [freddi]() {return freddi->Mdisk();}},
			{"Rhot", "Rsun", [freddi]() {return cmToSun(freddi->R()[freddi->last()]);}},
			{"Cirrout", "float", [freddi]() {return freddi->Cirr()[freddi->last()];}},
			{"H2R", "float", [freddi]() {return freddi->Height()[freddi->last()] / freddi->R()[freddi->last()];}},
			{"Teffout", "K", [freddi]() {return freddi->Tph()[freddi->last()];}},
			{"Tirrout", "K", [freddi]() {return freddi->Tirr()[freddi->last()];}},
			{"Qiir2Qvisout", "float", [freddi]() {return pow(freddi->Tirr()[freddi->last()] / freddi->Tph_vis()[freddi->last()], 4.);}},
			{"Lx", "erg/s", [freddi]() {return freddi->Lx();}},
			{"Lbol", "erg/s", [freddi]() {return freddi->Lbol_disk();}},
			{"Fx", "erg/s", [freddi]() {return freddi->Lx() * 2.0 * freddi->cosiOverD2() / FOUR_M_PI;}},
			{"mU", "mag", [freddi]() {return freddi->mU();}},
			{"mB", "mag", [freddi]() {return freddi->mB();}},
			{"mV", "mag", [freddi]() {return freddi->mV();}},
			{"mR", "mag", [freddi]() {return freddi->mR();}},
			{"mI", "mag", [freddi]() {return freddi->mI();}},
			{"mJ", "mag", [freddi]() {return freddi->mJ();}},
	};
	const bool cold_disk = freddi->args().flux->cold_disk;
	const bool star = freddi->args().flux->star;
	const auto& lambdas = freddi->args().flux->lambdas;
	for (size_t i = 0; i < lambdas.size(); ++i) {
		const double lambda = lambdas[i];
		fields.emplace_back(
				std::string("Fnu") + std::to_string(i),
				"erg/s/cm^2/Hz",
				[freddi, lambda]() { return freddi->flux(lambda); }
		);
		if (cold_disk) {
			fields.emplace_back(
					std::string("Fnu") + std::to_string(i) + "_cold",
					"erg/s/cm^2/Hz",
					[freddi, lambda]() { return freddi->flux_region<FreddiState::ColdRegion>(lambda); }
			);
		}
		if (star) {
			fields.emplace_back(
					"Fnu" + std::to_string(i) + "_star_min",
					"erg/s/cm^2/Hz",
					[freddi, lambda]() { return freddi->flux_star(lambda, 0.0); }
			);
			fields.emplace_back(
					"Fnu" + std::to_string(i) + "_star_max",
					"erg/s/cm^2/Hz",
					[freddi, lambda]() { return freddi->flux_star(lambda, M_PI); }
			);
		}
	}
	for (const auto& pb : freddi->args().flux->passbands) {
		fields.emplace_back(
				"Fnu" + pb.name,
				"erg/s/cm^2/Hz",
				[freddi, &pb]() { return freddi->flux(pb); }
		);
		if (cold_disk) {
			fields.emplace_back(
					std::string("Fnu") + pb.name + "_cold",
					"erg/s/cm^2/Hz",
					[freddi, &pb]() { return freddi->flux_region<FreddiState::ColdRegion>(pb); }
			);
		}
		if (star) {
			fields.emplace_back(
					"Fnu" + pb.name + "_star_min",
					"erg/s/cm^2/Hz",
					[freddi, &pb]() { return freddi->flux_star(pb, 0.0); }
			);
			fields.emplace_back(
					"Fnu" + pb.name + "_star_max",
					"erg/s/cm^2/Hz",
					[freddi, &pb]() { return freddi->flux_star(pb, M_PI); }
			);
		}
	}
	return fields;
}

std::vector<FileOutputLongField> FreddiFileOutput::initializeDiskStructureFields(const std::shared_ptr<FreddiEvolution>& freddi) {
	return {
			{"h", "cm^2/s", [freddi](size_t i) {return freddi->h()[i];}},
			{"R", "cm", [freddi](size_t i) {return freddi->R()[i];}},
			{"F", "dyn*cm", [freddi](size_t i) {return freddi->F()[i];}},
			{"Sigma", "g/cm^2", [freddi](size_t i) {return freddi->Sigma()[i];}},
			{"Teff", "K", [freddi](size_t i) {return freddi->Tph()[i];}},
			{"Tvis", "K", [freddi](size_t i) {return freddi->Tph_vis()[i];}},
			{"Tirr", "K", [freddi](size_t i) {return freddi->Tirr()[i];}},
			{"Height", "cm", [freddi](size_t i) {return freddi->Height()[i];}},
	};
}

std::vector<FileOutputLongField> FreddiFileOutput::initializeStarFields(const std::shared_ptr<FreddiEvolution>& freddi) {
	auto& star = freddi->star();
	const auto& triangles = star.triangles();

	std::vector<FileOutputLongField> fields;

	for (size_t i_vertex = 0; i_vertex < 3; ++i_vertex) {
		const auto prefix = "vertex" + std::to_string(i_vertex) + "_";
		fields.emplace_back(prefix + "x", "cm", [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].x();});
		fields.emplace_back(prefix + "y", "cm", [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].y();});
		fields.emplace_back(prefix + "z", "cm", [triangles, i_vertex](size_t i) {return triangles[i].vertices()[i_vertex].z();});
	}

	fields.emplace_back("center_x", "cm", [triangles](size_t i) {return triangles[i].center().x();});
	fields.emplace_back("center_y", "cm", [triangles](size_t i) {return triangles[i].center().y();});
	fields.emplace_back("center_z", "cm", [triangles](size_t i) {return triangles[i].center().z();});

	fields.emplace_back("Tth", "K", [freddi](size_t i) {return freddi->star().Tth()[i];});
	fields.emplace_back("Teff", "K", [freddi](size_t i) {return freddi->star().Teff()[i];});

	return fields;
}
