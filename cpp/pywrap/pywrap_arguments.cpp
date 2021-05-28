#include <string>

#include <boost/python/make_constructor.hpp>

#include <util.hpp>

#include "converters.hpp"
#include "pywrap_arguments.hpp"
#include "util.hpp"


using namespace boost::python;


boost::shared_ptr<GeneralArguments> make_general_arguments() {
	return boost::make_shared<GeneralArguments>("", "", 0, 1, false, false);
}

boost::shared_ptr<BasicDiskBinaryArguments> make_basic_disk_binary_arguments(
		double alpha, const object& alphacold,
		double Mx, double kerr,
		double period,
		double Mopt, double roche_lobe_fill, double Topt,
		const object& rin, const object& rout, const object& risco) {
	return boost::make_shared<BasicDiskBinaryArguments>(
			alpha, objToOpt<double>(alphacold),
			Mx, kerr,
			period,
			Mopt, roche_lobe_fill, Topt,
			objToOpt<double>(rin), objToOpt<double>(rout), objToOpt<double>(risco));
}

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot, double Tirr2Tvishot,
		const std::string& initialcond,
		const object& F0, const object& Mdisk0, const object& Mdot0,
		const object& powerorder, const object& gaussmu, const object& gausssigma,
		const std::string& wind, const object& windparams) {
	return boost::make_shared<DiskStructureArguments>(
			basic_disk_binary_arguments,
			opacity,
			Mdotout,
			boundcond, Thot, Tirr2Tvishot,
			initialcond,
			objToOpt<double>(F0), objToOpt<double>(Mdisk0), objToOpt<double>(Mdot0),
			objToOpt<double>(powerorder), objToOpt<double>(gaussmu), objToOpt<double>(gausssigma),
			wind, mapping_to_map(windparams));
}

boost::shared_ptr<SelfIrradiationArguments> make_self_irradiation_arguments(
		double Cirr, double irrindex,
		double Cirr_cold, double irrindex_cold, double height_to_radius_cold,
		const std::string& angular_dist_disk
		) {
	return boost::make_shared<SelfIrradiationArguments>(
			Cirr, irrindex,
			Cirr_cold, irrindex_cold, height_to_radius_cold,
			angular_dist_disk);
}

boost::shared_ptr<FluxArguments> make_flux_arguments(
		double colourfactor,
		double emin, double emax,
		double star_albedo,
		double inclination, double ephemeris_t0, double distance) {
	return boost::make_shared<FluxArguments>(
			colourfactor,
			emin, emax,
			star_albedo,
			inclination, ephemeris_t0, distance,
			false, false,
			vecd(),
			std::vector<Passband>());
}

boost::shared_ptr<CalculationArguments> make_calculation_arguments(
		double inittime,
		double time, const object& tau,
		unsigned int Nx, const std::string& gridscale, const unsigned short starlod,
		const object& eps) {
	if (eps.ptr() == object().ptr()) {
		return boost::make_shared<CalculationArguments>(inittime, time, objToOpt<double>(tau), Nx, gridscale, starlod);
	}
	return boost::make_shared<CalculationArguments>(inittime, time, objToOpt<double>(tau), Nx, gridscale, starlod, extract<double>(eps));
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
		const std::string& fptype, const object& fpparams_,
		const std::string& kappat_type, const object& kappat_params_,
		const std::string& ns_grav_redshift) {
	return boost::make_shared<NeutronStarArguments>(
			nsprop,
			objToOpt<double>(freqx), objToOpt<double>(Rx), Bx, hotspotarea,
			epsilonAlfven, inversebeta, Rdead,
			fptype, mapping_to_map(fpparams_),
			kappat_type, mapping_to_map(kappat_params_),
			ns_grav_redshift);
}

boost::shared_ptr<NeutronStarBasicDiskBinaryArguments> make_neutron_star_basic_disk_binary_arguments(
		const NeutronStarArguments& ns_args,
		double alpha, const object& alphacold,
		double Mx, double kerr,
		double period,
		double Mopt, double roche_lobe_fill, double Topt,
		const object& rin, const object& rout, const object& risco) {
	return boost::make_shared<NeutronStarBasicDiskBinaryArguments>(
			ns_args,
			alpha, objToOpt<double>(alphacold),
			Mx, kerr,
			period,
			Mopt, roche_lobe_fill, Topt,
			objToOpt<double>(rin), objToOpt<double>(rout), objToOpt<double>(risco));
}

boost::shared_ptr<NeutronStarSelfIrradiationArguments> make_neutron_star_self_irradiation_arguments(
		double Cirr, double irrindex,
		double Cirr_cold, double irrindex_cold, double height_to_radius_cold,
		const std::string& angular_dist_disk, const std::string& angular_dist_ns) {
	return boost::make_shared<NeutronStarSelfIrradiationArguments>(
			Cirr, irrindex,
			Cirr_cold, irrindex_cold, height_to_radius_cold,
			angular_dist_disk, angular_dist_ns);
}

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const GeneralArguments& general,
		const NeutronStarBasicDiskBinaryArguments& basic,
		const DiskStructureArguments& disk,
		const NeutronStarSelfIrradiationArguments& irr,
		const FluxArguments& flux,
		const CalculationArguments& calc,
		const NeutronStarArguments& ns) {
	return boost::make_shared<FreddiNeutronStarArguments>(
			new GeneralArguments(general),
			new NeutronStarBasicDiskBinaryArguments(basic),
			new DiskStructureArguments(disk),
			new NeutronStarSelfIrradiationArguments(irr),
			new FluxArguments(flux),
			new CalculationArguments(calc),
			new NeutronStarArguments(ns));
}


void wrap_arguments() {
	class_<GeneralArguments>("_GeneralArguments", no_init)
		.def("__init__", make_constructor(make_general_arguments))
	;
	auto _BasicDiskBinaryArguments = class_<BasicDiskBinaryArguments>("_BasicDiskBinaryArguments", no_init)
		.def("__init__", make_constructor(make_basic_disk_binary_arguments))
	;
	class_<DiskStructureArguments>("_DiskStructureArguments", no_init)
		.def("__init__", make_constructor(make_disk_structure_arguments))
	;
	class_<SelfIrradiationArguments>("_SelfIrradiationArguments", no_init)
		.def("__init__", make_constructor(make_self_irradiation_arguments))
	;
	class_<FluxArguments>("_FluxArguments", no_init)
		.def("__init__", make_constructor(make_flux_arguments))
	;
	class_<CalculationArguments>("_CalculationArguments", no_init)
		.def("__init__", make_constructor(make_calculation_arguments))
	;
	class_<FreddiArguments>("_Arguments", no_init)
		.def("__init__", make_constructor(make_freddi_arguments))
	;
	class_<NeutronStarArguments>("_NeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_neutron_star_arguments))
	;
	class_< FreddiNeutronStarArguments, bases<FreddiArguments> >("_FreddiNeutronStarArguments", no_init)
		.def("__init__", make_constructor(make_freddi_neutron_star_arguments))
	;
}
