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
	return fields;
}

std::vector<FileOutputLongField> FreddiNeutronStarFileOutput::initializeLongFields(const std::shared_ptr<FreddiNeutronStarEvolution>& freddi) {
	auto fields = FreddiFileOutput::initializeLongFields(freddi);
	fields.emplace_back("Fmagn", "dyn*cm", [freddi]() {return freddi->Fmagn();});
	return fields;
}
