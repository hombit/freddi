#include <iostream>

#include <ns/ns_output.hpp>


std::vector<FileOutputShortField> FreddiNeutronStarFileOutput::initializeShortFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeShortFields(freddi);
	fields.emplace_back("Rin", "cm", [freddi]() {return freddi->R()[freddi->first()];});
	fields.emplace_back("Lxns", "erg/s", [freddi]() {return freddi->Lx_ns();});
	fields.emplace_back("Lbolns", "erg/s", [freddi]() {return freddi->Lbol_ns();});
	fields.emplace_back("Thotspot", "keV", [freddi]() {return kToKev(freddi->T_hot_spot());});
	fields.emplace_back("fpin", "float", [freddi]() {return freddi->fp();});
	fields.emplace_back("Fmagnin", "dyn*cm", [freddi]() {return freddi->Fmagn()[freddi->first()];});
	fields.emplace_back("Fin", "dyn*cm", [freddi]() {return freddi->F()[freddi->first()];});
	for (const auto& pb : freddi->args().flux->passbands) {
		fields.emplace_back(
				std::string("Fnu") + pb.name + "_star",
				"erg/s/cm^2/Hz",
				[freddi, &pb]() {
					const double Tth = 3500.0;
					const auto& basic = freddi->args().basic;
					const double Rstar = BinaryFunctions::rocheLobeVolumeRadius(basic->Mopt, basic->Mx, basic->period);
					const double semiaxis = BinaryFunctions::semiaxis(basic->Mx, basic->Mopt, basic->period);
					const double Rstar2semiaxis = BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(basic->Mopt / basic->Mx);
					const double star_area = M_PI * m::pow<2>(Rstar);
					const double Lbol_disk = (freddi->F()[freddi->first()] + 0.5 * freddi->Mdot_in() * freddi->h()[freddi->first()]) * freddi->omega_i(freddi->first());
					const double Qirr = (freddi->Lbol_ns() + Lbol_disk * Rstar2semiaxis) / (4.0 * M_PI * m::pow<2>(semiaxis));
					std::cout << Rstar << std::endl;
					std::cout << Qirr << " " << freddi->Qx()[freddi->Nx()-1] << std::endl;
					const double Tstarspot = std::pow(Qirr / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT + m::pow<4>(Tth), 0.25);
//					const double F = star_area * (pb.bb_nu(Tstarspot) - pb.bb_nu(Tth)) / m::pow<2>(freddi->distance());
					const double F = star_area * pb.bb_nu(Tth) / m::pow<2>(freddi->distance());
					return F;
				}
		);
	}
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeLongFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeLongFields(freddi);
	fields.emplace_back("Fmagn", "dyn*cm", [freddi]() {return freddi->Fmagn();});
	return fields;
}
