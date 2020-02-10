#include <string>

#include <boost/python/make_constructor.hpp>

#include <util.hpp>

#include "converters.hpp"
#include "pywrap_arguments.hpp"
#include "util.hpp"


using namespace boost::python;


boost::shared_ptr<GeneralArguments> make_general_arguments() {
	return boost::make_shared<GeneralArguments>("", "", false);
}

boost::shared_ptr<BasicDiskBinaryArguments> make_basic_disk_binary_arguments(
		double alpha,
		double Mx, double kerr,
		double period,
		double Mopt, const object& Ropt, double Topt,
		const object& rin, const object& rout, const object& risco) {
	return boost::make_shared<BasicDiskBinaryArguments>(
			alpha,
			Mx, kerr,
			period,
			Mopt, objToOpt<double>(Ropt), Topt,
			objToOpt<double>(rin), objToOpt<double>(rout), objToOpt<double>(risco));
}

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot,
		const std::string& initialcond,
		const object& F0, const object& Mdisk0, const object& Mdot0,
		const object& powerorder, const object& gaussmu, const object& gausssigma,
		const std::string& wind, const object& windparams) {
	return boost::make_shared<DiskStructureArguments>(
			basic_disk_binary_arguments,
			opacity,
			Mdotout,
			boundcond, Thot,
			initialcond,
			objToOpt<double>(F0), objToOpt<double>(Mdisk0), objToOpt<double>(Mdot0),
			objToOpt<double>(powerorder), objToOpt<double>(gaussmu), objToOpt<double>(gausssigma),
			wind, mapping_to_map(windparams));
}

boost::shared_ptr<SelfIrradiationArguments> make_self_irradiation_arguments(double Cirr, const std::string& irrfactortype) {
	return boost::make_shared<SelfIrradiationArguments>(Cirr, irrfactortype);
}

boost::shared_ptr<FluxArguments> make_flux_arguments(
		double colourfactor,
		double emin, double emax,
		double star_albedo,
		double inclination, double distance) {
	return boost::make_shared<FluxArguments>(
			colourfactor,
			emin, emax,
			star_albedo,
			inclination, distance,
			false, false,
			vecd(),
			std::vector<Passband>());
}

boost::shared_ptr<CalculationArguments> make_calculation_arguments(
		double time, const object& tau,
		unsigned int Nx, const std::string& gridscale, const unsigned short starlod,
		const object& eps) {
	if (eps.ptr() == object().ptr()) {
		return boost::make_shared<CalculationArguments>(time, objToOpt<double>(tau), Nx, gridscale, starlod);
	}
	return boost::make_shared<CalculationArguments>(time, objToOpt<double>(tau), Nx, gridscale, starlod, extract<double>(eps));
}

boost::shared_ptr<FreddiArguments> make_freddi_arguments(
		const GeneralArguments& general,
		const BasicDiskBinaryArguments& basic,
		const DiskStructureArguments& disk,
		const SelfIrradiationArguments& irr,
		const FluxArguments& flux,
		const CalculationArguments& calc) {
	return boost::make_shared<FreddiArguments>(
			new GeneralArguments(general),
			new BasicDiskBinaryArguments(basic),
			new DiskStructureArguments(disk),
			new SelfIrradiationArguments(irr),
			new FluxArguments(flux),
			new CalculationArguments(calc));
}

boost::shared_ptr<NeutronStarArguments> make_neutron_star_arguments(
		const std::string& nsprop,
		const object& freqx, const object& Rx, double Bx, double hotspotarea,
		double epsilonAlfven, double inversebeta, double Rdead,
		const std::string& fptype, const object& fpparams_) {
	return boost::make_shared<NeutronStarArguments>(
			nsprop,
			objToOpt<double>(freqx), objToOpt<double>(Rx), Bx, hotspotarea,
			epsilonAlfven, inversebeta, Rdead,
			fptype, mapping_to_map(fpparams_));
}

boost::shared_ptr<NeutronStarBasicDiskBinaryArguments> make_neutron_star_basic_disk_binary_arguments(
		const NeutronStarArguments& ns_args,
		double alpha,
		double Mx, double kerr,
		double period,
		double Mopt, const object& Ropt, double Topt,
		const object& rin, const object& rout, const object& risco) {
	return boost::make_shared<NeutronStarBasicDiskBinaryArguments>(
			ns_args,
			alpha,
			Mx, kerr,
			period,
			Mopt, objToOpt<double>(Ropt), Topt,
			objToOpt<double>(rin), objToOpt<double>(rout), objToOpt<double>(risco));
}

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const GeneralArguments& general,
		const NeutronStarBasicDiskBinaryArguments& basic,
		const DiskStructureArguments& disk,
		const SelfIrradiationArguments& irr,
		const FluxArguments& flux,
		const CalculationArguments& calc,
		const NeutronStarArguments& ns) {
	return boost::make_shared<FreddiNeutronStarArguments>(
			new GeneralArguments(general),
			new NeutronStarBasicDiskBinaryArguments(basic),
			new DiskStructureArguments(disk),
			new SelfIrradiationArguments(irr),
			new FluxArguments(flux),
			new CalculationArguments(calc),
			new NeutronStarArguments(ns));
}


void wrap_arguments() {
	class_<GeneralArguments>("_GeneralArguments", no_init)
		.def("__init__", make_constructor(make_general_arguments))
	;
	auto _BasicDiskBinaryArguments = class_<BasicDiskBinaryArguments>("_BasicDiskBinaryArguments", no_init)
		.def("__init__", make_constructor(make_basic_disk_binary_arguments, default_call_policies(),
				(arg("alpha"),
				 arg("Mx"),
				 arg("kerr"),
				 arg("Mopt"),
				 arg("Topt"),
				 arg("period"),
				 arg("rin"),
				 arg("rout")))
			)
	;
	class_<DiskStructureArguments>("_DiskStructureArguments", no_init)
		.def("__init__", make_constructor(make_disk_structure_arguments, default_call_policies(),
				(arg("basic_disk_binary_arguments"),
				 arg("opacity"),
				 arg("Mdotout"),
				 arg("boundcond"),
				 arg("Thot"),
				 arg("initialcond"),
				 arg("F0"),
				 arg("powerorder"),
				 arg("gaussmu"),
				 arg("gausssigma"),
				 arg("Mdisk0"),
				 arg("Mdot0"),
				 arg("wind"),
				 arg("windparams")))
			)
	;
	class_<SelfIrradiationArguments>("_SelfIrradiationArguments", no_init)
		.def("__init__", make_constructor(make_self_irradiation_arguments, default_call_policies(),
				(arg("Cirr"),
				 arg("irrfactortype")))
			)
	;
	class_<FluxArguments>("_FluxArguments", no_init)
		.def("__init__", make_constructor(make_flux_arguments, default_call_policies(),
				(arg("colourfactor"),
				 arg("emin"),
				 arg("emax"),
				 arg("inclination"),
				 arg("distance")))
			)
	;
	class_<CalculationArguments>("_CalculationArguments", no_init)
		.def("__init__", make_constructor(make_calculation_arguments, default_call_policies(),
				(arg("time"),
				 arg("tau"),
				 arg("Nx"),
				 arg("gridscale"),
				 arg("starlod"),
				 arg("eps")))
			)
	;
	class_<FreddiArguments>("_Arguments", no_init)
		.def("__init__", make_constructor(make_freddi_arguments))
	;
	class_<NeutronStarArguments>("_NeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_neutron_star_arguments, default_call_policies(),
				(arg("nsprop"),
				 arg("Rx"),
				 arg("freqx"),
				 arg("Bx"),
				 arg("hotspotarea"),
				 arg("epsilonAlfven"),
				 arg("inversebeta"),
				 arg("Rdead"),
				 arg("fptype"),
				 arg("fpparams")))
			)
	;
	class_< FreddiNeutronStarArguments, bases<FreddiArguments> >("_FreddiNeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_freddi_neutron_star_arguments))
	;
}
