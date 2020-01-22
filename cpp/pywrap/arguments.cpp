#include <string>

#include <boost/python/make_constructor.hpp>

#include <util.hpp>

#include "arguments.hpp"
#include "converters.hpp"


using namespace boost::python;


boost::shared_ptr<GeneralArguments> make_general_arguments() {
	return boost::make_shared<GeneralArguments>("", "", false);
}

boost::shared_ptr<BasicDiskBinaryArguments> make_basic_disk_binary_arguments(
		double alpha,
		double Mx, double kerr,
		double Mopt, double Topt,
		double period,
		const object& rin, const object& rout) {
	const auto None = object().ptr();
	if (rin.ptr() == None && rout.ptr() == None) {
		return boost::make_shared<BasicDiskBinaryArguments>(BasicDiskBinaryArguments::constructWithoutRinRout(alpha, Mx, kerr, Mopt, Topt, period));
	}
	if (rin.ptr() == None) {
		return boost::make_shared<BasicDiskBinaryArguments>(BasicDiskBinaryArguments::constructWithoutRin(alpha, Mx, kerr, Mopt, Topt, period, extract<double>(rout)));
	}
	if (rout.ptr() == None) {
		return boost::make_shared<BasicDiskBinaryArguments>(BasicDiskBinaryArguments::constructWithoutRout(alpha, Mx, kerr, Mopt, Topt, period, extract<double>(rin)));
	}
	return boost::make_shared<BasicDiskBinaryArguments>(alpha, Mx, kerr, Mopt, Topt, period, extract<double>(rin), extract<double>(rout));
}

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot,
		const std::string& initialcond, double F0,
		double powerorder, double gaussmu, double gausssigma,
		const object& Mdisk0_, const object& Mdot0_,
		const std::string& wind, const object& windparams_) {
	const auto None = object().ptr();

	double Mdisk0 = -1;
	bool is_Mdisk0_specified = false;
	if (Mdisk0_.ptr() != None) {
		Mdisk0 = extract<double>(Mdisk0_);
		is_Mdisk0_specified = true;
	}

	double Mdot0 = -1;
	bool is_Mdot0_specified = false;
	if (Mdot0_.ptr() != None) {
		Mdot0 = extract<double>(Mdot0_);
		is_Mdot0_specified = true;
	}

	return boost::make_shared<DiskStructureArguments>(
			basic_disk_binary_arguments,
			opacity,
			Mdotout,
			boundcond, Thot,
			initialcond, F0,
			powerorder, gaussmu, gausssigma,
			is_Mdisk0_specified, is_Mdot0_specified,
			Mdisk0, Mdot0,
			wind, mapping_to_map(windparams_));
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
		double time, double tau, unsigned int Nx, const std::string& gridscale, const unsigned short starlod,
		const object& eps) {
	if (eps.ptr() == object().ptr()) {
		return boost::make_shared<CalculationArguments>(time, tau, Nx, gridscale, starlod);
	}
	return boost::make_shared<CalculationArguments>(time, tau, Nx, gridscale, starlod, extract<double>(eps));
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
		double Rx, double freqx, double Bx, double hotspotarea,
		double epsilonAlfven, double inversebeta, double Rdead,
		const std::string& fptype, const object& fpparams_) {
	return boost::make_shared<NeutronStarArguments>(
			Rx, freqx, Bx, hotspotarea,
			epsilonAlfven, inversebeta, Rdead,
			fptype, mapping_to_map(fpparams_));
}

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const FreddiArguments& freddi_args, const NeutronStarArguments& ns) {
	return boost::make_shared<FreddiNeutronStarArguments>(freddi_args, new NeutronStarArguments(ns));
}


void wrap_arguments() {
	class_<GeneralArguments>("_GeneralArguments", no_init)
		.def("__init__", make_constructor(make_general_arguments))
	;
	auto _BasicDiskBinaryArguments = class_<BasicDiskBinaryArguments>("_BasicDiskBinaryArguments", no_init)
		.def("__init__", make_constructor(make_basic_disk_binary_arguments, default_call_policies(),
				(arg("alpha")=BasicDiskBinaryArguments::default_alpha,
				 arg("Mx")=BasicDiskBinaryArguments::default_Mx,
				 arg("kerr")=BasicDiskBinaryArguments::default_kerr,
				 arg("Mopt")=BasicDiskBinaryArguments::default_Mopt,
				 arg("Topt")=BasicDiskBinaryArguments::default_Topt,
				 arg("period")=BasicDiskBinaryArguments::default_period,
				 arg("rin")=object(),
				 arg("rout")=object()))
			)
		.def_readonly("default_alpha", &BasicDiskBinaryArguments::default_alpha)
	;
	class_<DiskStructureArguments>("_DiskStructureArguments", no_init)
		.def("__init__", make_constructor(make_disk_structure_arguments, default_call_policies(),
				(arg("basic_disk_binary_arguments"),
				 arg("opacity")=DiskStructureArguments::default_opacity,
				 arg("Mdotout")=DiskStructureArguments::default_Mdotout,
				 arg("boundcond")=DiskStructureArguments::default_boundcond,
				 arg("Thot")=DiskStructureArguments::default_Thot,
				 arg("initialcond")=DiskStructureArguments::default_initialcond,
				 arg("F0")=DiskStructureArguments::default_F0,
				 arg("powerorder")=DiskStructureArguments::default_powerorder,
				 arg("gaussmu")=DiskStructureArguments::default_gaussmu,
				 arg("gausssigma")=DiskStructureArguments::default_gausssigma,
				 arg("Mdisk0")=object(),
				 arg("Mdot0")=object(),
				 arg("wind")=DiskStructureArguments::default_wind,
				 arg("windparams")=tuple()))
			)
	;
	class_<SelfIrradiationArguments>("_SelfIrradiationArguments", no_init)
		.def("__init__", make_constructor(make_self_irradiation_arguments, default_call_policies(),
				(arg("Cirr")=SelfIrradiationArguments::default_Cirr,
				 arg("irrfactortype")=SelfIrradiationArguments::default_irrfactortype))
			)
	;
	class_<FluxArguments>("_FluxArguments", no_init)
		.def("__init__", make_constructor(make_flux_arguments, default_call_policies(),
				(arg("colourfactor")=FluxArguments::default_colourfactor,
				 arg("emin")=FluxArguments::default_emin,
				 arg("emax")=FluxArguments::default_emax,
				 arg("inclination")=FluxArguments::default_inclination,
				 arg("distance")=FluxArguments::default_distance))
			)
	;
	class_<CalculationArguments>("_CalculationArguments", no_init)
		.def("__init__", make_constructor(make_calculation_arguments, default_call_policies(),
				(arg("time")=CalculationArguments::default_time,
				 arg("tau")=CalculationArguments::default_tau,
				 arg("Nx")=CalculationArguments::default_Nx,
				 arg("gridscale")=CalculationArguments::default_gridscale,
				 arg("starlod")=CalculationArguments::default_starlod,
				 arg("eps")=object()))
			)
	;
	class_<FreddiArguments>("_Arguments", no_init)
		.def("__init__", make_constructor(make_freddi_arguments))
	;
	class_<NeutronStarArguments>("_NeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_neutron_star_arguments, default_call_policies(),
				(arg("Rx")=NeutronStarArguments::default_Rx,
				 arg("freqx")=NeutronStarArguments::default_freqx,
				 arg("Bx")=NeutronStarArguments::default_Bx,
				 arg("hotspotarea")=NeutronStarArguments::default_hotspotarea,
				 arg("epsilonAlfven")=NeutronStarArguments::default_epsilonAlfven,
				 arg("inversebeta")=NeutronStarArguments::default_inversebeta,
				 arg("Rdead")=NeutronStarArguments::default_Rdead,
				 arg("fptype")=NeutronStarArguments::default_fptype,
				 arg("fpparams")=tuple()))
			)
	;
	class_< FreddiNeutronStarArguments, bases<FreddiArguments> >("_FreddiNeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_freddi_neutron_star_arguments))
	;
}
