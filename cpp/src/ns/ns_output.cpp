#include <iostream>

#include <ns/ns_output.hpp>


std::vector<FileOutputShortField> FreddiNeutronStarFileOutput::initializeShortFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeShortFields(freddi);
	fields.emplace_back("Rin", "cm", "Inner radius of the disk", [freddi]() {return freddi->R()[freddi->first()];});
	fields.emplace_back("etans", "float", "Accretion efficiency of the neutron star", [freddi]() {return freddi->eta_ns();});
	fields.emplace_back("Lxns", "erg/s", "X-ray luminosity of the neutron star in the given energy range [emin, emax]", [freddi]() {return freddi->Lx_ns();});
	fields.emplace_back("Lbolns", "erg/s", "Bolometric luminosity of the neutron star", [freddi]() {return freddi->Lbol_ns();});
	fields.emplace_back("Fxns", "erg/s/cm^2", "X-ray flux of the neutron star in the given energy range [emin, emax]", [freddi]() {return freddi->Lx_ns() * freddi->angular_dist_ns(freddi->cosi()) / (FOUR_M_PI * m::pow<2>(freddi->distance()));}),
	fields.emplace_back("Fbolns", "erg/s/cm^2", "Bolometric flux of the neutron star", [freddi]() {return freddi->Lbol_ns() * freddi->angular_dist_ns(freddi->cosi()) / (FOUR_M_PI * m::pow<2>(freddi->distance()));}),
	fields.emplace_back("Thotspot", "keV", "Temperature of the neutron star 'hot spot'", [freddi]() {return kToKev(freddi->T_hot_spot());});
	fields.emplace_back("fpin", "float", "Part of accreting matter falling onto the neutron star", [freddi]() {return freddi->fp();});
	fields.emplace_back("Fmagnin", "dyn*cm", "Magnetic torque at the inner radius of the disk", [freddi]() {return freddi->Fmagn()[freddi->first()];});
	fields.emplace_back("Fin", "dyn*cm", "Viscous torque at the inner radius of the disk", [freddi]() {return freddi->F()[freddi->first()];});
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeDiskStructureFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeDiskStructureFields(freddi);
	fields.emplace_back("Fmagn", "dyn*cm", "Magnetic torque", [freddi](size_t i) {return freddi->Fmagn()[i];});
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeStarFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	return FreddiFileOutput::initializeStarFields(freddi);
}
