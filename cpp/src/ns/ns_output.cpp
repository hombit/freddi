#include <iostream>

#include <ns/ns_output.hpp>


std::vector<FileOutputShortField> FreddiNeutronStarFileOutput::initializeShortFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeShortFields(freddi);
	fields.emplace_back("Rin", "cm", [freddi]() {return freddi->R()[freddi->first()];});
	fields.emplace_back("Lxns", "erg/s", [freddi]() {return freddi->Lx_ns();});
	fields.emplace_back("Lbolns", "erg/s", [freddi]() {return freddi->Lbol_ns();});
	fields.emplace_back("Fxns", "erg/s", [freddi]() {return freddi->Lx_ns() / (FOUR_M_PI * m::pow<2>(freddi->distance()));}),
	fields.emplace_back("Thotspot", "keV", [freddi]() {return kToKev(freddi->T_hot_spot());});
	fields.emplace_back("fpin", "float", [freddi]() {return freddi->fp();});
	fields.emplace_back("Fmagnin", "dyn*cm", [freddi]() {return freddi->Fmagn()[freddi->first()];});
	fields.emplace_back("Fin", "dyn*cm", [freddi]() {return freddi->F()[freddi->first()];});
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeDiskStructureFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeDiskStructureFields(freddi);
	fields.emplace_back("Fmagn", "dyn*cm", [freddi](size_t i) {return freddi->Fmagn()[i];});
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeStarFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	return FreddiFileOutput::initializeStarFields(freddi);
}
