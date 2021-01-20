#include <sstream>

#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>

#include <freddi_evolution.hpp>
#include <unit_transformation.hpp>
#include <ns/ns_evolution.hpp>

#include "pywrap_arguments.hpp"
#include "pywrap_freddi_evolution.hpp"

using namespace boost::python;

dict evolution_required_args() {
	dict kw;

	kw["alpha"] = 0.0;
	kw["Mx"] = 0.0;
	kw["Mopt"] = 0.0;
	kw["period"] = 0.0;

	kw["initialcond"] = "quasistat";
	kw["F0"] = 0.0;

	kw["distance"] = 0.0;

	kw["time"] = 0.0;

	return kw;
}

dict evolution_kwdefaults() {
	dict kw;

	kw["__cgs"] = true;

	kw["alphacold"] = object();
	kw["kerr"] = BasicDiskBinaryArguments::default_kerr;
	kw["rochelobefill"] = BasicDiskBinaryArguments::default_roche_lobe_fill;
	kw["Topt"] = BasicDiskBinaryArguments::default_Topt;
	kw["rin"] = object();
	kw["rout"] = object();
	kw["risco"] = object();

	kw["opacity"] = DiskStructureArguments::default_opacity;
	kw["Mdotout"] = DiskStructureArguments::default_Mdotout;
	kw["boundcond"] = DiskStructureArguments::default_boundcond;
	kw["initialcond"] = DiskStructureArguments::default_initialcond;
	kw["Thot"] = DiskStructureArguments::default_Thot;
	kw["Qirr2Qvishot"] = m::pow<4>(DiskStructureArguments::default_Tirr2Tvishot);
	kw["F0"] = object();
	kw["Mdisk0"] = object();
	kw["Mdot0"] = object();
	kw["powerorder"] = object();
	kw["gaussmu"] = object();
	kw["gausssigma"] = object();
	kw["wind"] = DiskStructureArguments::default_wind;
	kw["windparams"] = tuple();

	kw["Cirr"] = SelfIrradiationArguments::default_Cirr;
	kw["irrindex"] = SelfIrradiationArguments::default_irrindex;
	kw["Cirrcold"] = SelfIrradiationArguments::default_Cirr_cold;
	kw["irrindexcold"] = SelfIrradiationArguments::default_irrindex_cold;
	kw["h2rcold"] = SelfIrradiationArguments::default_height_to_radius_cold;
	kw["angulardistdisk"] = SelfIrradiationArguments::default_angular_dist_disk;

	kw["colourfactor"] = FluxArguments::default_colourfactor;
	kw["emin"] = FluxArguments::default_emin;
	kw["emax"] = FluxArguments::default_emax;
	kw["staralbedo"] = FluxArguments::default_star_albedo;
	kw["inclination"] = FluxArguments::default_inclination;
	kw["ephemerist0"] = FluxArguments::default_ephemeris_t0;

	kw["inittime"] = CalculationArguments::default_init_time;
	kw["tau"] = object();
	kw["Nx"] = CalculationArguments::default_Nx;
	kw["gridscale"] = CalculationArguments::default_gridscale;
	kw["starlod"] = CalculationArguments::default_starlod;
	kw["eps"] = object();

	return kw;
}


void check_args(const tuple &args) {
	if (len(args) > 1) {
		PyErr_SetString(PyExc_TypeError, "Positional arguments are denied");
		throw error_already_set();
	}
}


void check_kwargs(const tuple& args, const dict& kwargs, const dict& required_args, const dict& kwdefaults) {
	stl_input_iterator<object> begin(kwargs), end;
	for (auto key = begin; key != end; ++key) {
		if (! kwdefaults.has_key(*key) && ! required_args.has_key(*key)) {
			std::stringstream msg;
			msg << extract<char*>(str(args[0])) << " got an unexpected keyword argument " << extract<char*>(*key);
			PyErr_SetString(PyExc_TypeError, msg.str().data());
			throw error_already_set();
		}
	}
}


void update_no_cgs_kwargs(dict& kw) {
	if (extract<bool>(kw["__cgs"])) return;

	static const auto None = object().ptr();

	kw["Mx"] = sunToGram(extract<double>(kw["Mx"]));
	kw["Mopt"] = sunToGram(extract<double>(kw["Mopt"]));
	kw["period"] = dayToS(extract<double>(kw["period"]));
	if (object(kw["rin"]).ptr() != None) {
		kw["rin"] = rgToCm(extract<double>(kw["rin"]), extract<double>(kw["Mx"]));
	}
	if (object(kw["rout"]).ptr() != None) {
		kw["rout"] = sunToCm(extract<double>(kw["rout"]));
	}
	kw["emin"] = kevToHertz(extract<double>(kw["emin"]));
	kw["emax"] = kevToHertz(extract<double>(kw["emax"]));
	kw["distance"] = kpcToCm(extract<double>(kw["distance"]));

	kw["inittime"] = dayToS(extract<double>(kw["inittime"]));
	kw["time"] = dayToS(extract<double>(kw["time"]));
	if (object(kw["tau"]).ptr() != None) {
		kw["tau"] = dayToS(extract<double>(kw["tau"]));
	}
}


boost::shared_ptr<FreddiArguments> make_freddi_arguments(dict& kw) {
	update_no_cgs_kwargs(kw);

	const auto general = make_general_arguments();
	const auto basic = make_basic_disk_binary_arguments(
			extract<double>(kw["alpha"]), kw["alphacold"],
			extract<double>(kw["Mx"]), extract<double>(kw["kerr"]),
			extract<double>(kw["period"]),
			extract<double>(kw["Mopt"]), extract<double>(kw["rochelobefill"]), extract<double>(kw["Topt"]),
			kw["rin"], kw["rout"], kw["risco"]);
	const auto disk = make_disk_structure_arguments(
			*basic,
			extract<std::string>(kw["opacity"]),
			extract<double>(kw["Mdotout"]),
			extract<std::string>(kw["boundcond"]), extract<double>(kw["Thot"]),
			std::pow(extract<double>(kw["Qirr2Qvishot"]), 0.25),
			extract<std::string>(kw["initialcond"]),
			kw["F0"], kw["Mdisk0"], kw["Mdot0"],
			kw["powerorder"], kw["gaussmu"], kw["gausssigma"],
			extract<std::string>(kw["wind"]), kw["windparams"]);
	const auto irr = make_self_irradiation_arguments(
			extract<double>(kw["Cirr"]),
			extract<double>(kw["irrindex"]),
			extract<double>(kw["Cirrcold"]),
			extract<double>(kw["irrindexcold"]),
			extract<double>(kw["h2rcold"]),
			extract<std::string>(kw["angulardistdisk"]));
	const auto flux = make_flux_arguments(
			extract<double>(kw["colourfactor"]),
			extract<double>(kw["emin"]), extract<double>(kw["emax"]),
			extract<double>(kw["staralbedo"]),
			extract<double>(kw["inclination"]), extract<double>(kw["ephemerist0"]), extract<double>(kw["distance"]));
	const auto calc = make_calculation_arguments(
			extract<double>(kw["inittime"]),
			extract<double>(kw["time"]), kw["tau"],
			extract<unsigned int>(kw["Nx"]), extract<std::string>(kw["gridscale"]),
			extract<unsigned short>(kw["starlod"]),
			kw["eps"]);
	// To avoid copy, create FreddiArguments constructor that accepts shared_ptr
	return make_freddi_arguments(*general, *basic, *disk, *irr, *flux, *calc);
}


object raw_make_evolution(const tuple& args, const dict& kwargs) {
	check_args(args);

	// It should be static, but segfault
	const auto kwdefaults = evolution_kwdefaults();
	check_kwargs(args, kwargs, evolution_required_args(), kwdefaults);

	dict kw;
	kw.update(kwdefaults);
	kw.update(kwargs);

	const auto freddi_args = make_freddi_arguments(kw);

	args[0].attr("_kwargs") = kw;
	return args[0].attr("__init__")(*freddi_args);
}


dict neutron_star_evolution_required_args() {
	auto kw = evolution_required_args();

	kw["Bx"] = 0.0;

	return kw;
}


dict neutron_star_evolution_kwdefaults() {
	dict kw;
	kw.update(evolution_kwdefaults());

	kw["angulardistns"] = NeutronStarSelfIrradiationArguments::default_angular_dist_ns;

	kw["nsprop"] = NeutronStarArguments::default_nsprop;
	kw["freqx"] = object();
	kw["Rx"] = object();
	kw["hotspotarea"] = NeutronStarArguments::default_hotspotarea;
	kw["epsilonAlfven"] = NeutronStarArguments::default_epsilonAlfven;
	kw["inversebeta"] = NeutronStarArguments::default_inversebeta;
	kw["Rdead"] = NeutronStarArguments::default_Rdead;
	kw["fptype"] = NeutronStarArguments::default_fptype;
	kw["fpparams"] = tuple();
	kw["kappattype"] = NeutronStarArguments::default_kappat_type;
	kw["kappatparams"] = dict();
	kw["kappatparams"]["value"] = kw["kappatparams"]["in"] = kw["kappatparams"]["out"] = NeutronStarArguments::default_kappat_value;
	kw["nsgravredshift"] = NeutronStarArguments::default_ns_grav_redshift;

	return kw;
}

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(dict& kw) {
	update_no_cgs_kwargs(kw);

	const auto ns_args = make_neutron_star_arguments(
			extract<std::string>(kw["nsprop"]),
			kw["freqx"], kw["Rx"], extract<double>(kw["Bx"]), extract<double>(kw["hotspotarea"]),
			extract<double>(kw["epsilonAlfven"]), extract<double>(kw["inversebeta"]), extract<double>(kw["Rdead"]),
			extract<std::string>(kw["fptype"]), kw["fpparams"],
			extract<std::string>(kw["kappattype"]), kw["kappatparams"],
			extract<std::string>(kw["nsgravredshift"]));

	const auto general = make_general_arguments();
	const auto basic = make_neutron_star_basic_disk_binary_arguments(
			*ns_args,
			extract<double>(kw["alpha"]), kw["alphacold"],
			extract<double>(kw["Mx"]), extract<double>(kw["kerr"]),
			extract<double>(kw["period"]),
			extract<double>(kw["Mopt"]), extract<double>(kw["rochelobefill"]), extract<double>(kw["Topt"]),
			kw["rin"], kw["rout"], kw["risco"]);
	const auto disk = make_disk_structure_arguments(
			*basic,
			extract<std::string>(kw["opacity"]),
			extract<double>(kw["Mdotout"]),
			extract<std::string>(kw["boundcond"]), extract<double>(kw["Thot"]),
			std::pow(extract<double>(kw["Qirr2Qvishot"]), 0.25),
			extract<std::string>(kw["initialcond"]),
			kw["F0"], kw["Mdisk0"], kw["Mdot0"],
			kw["powerorder"], kw["gaussmu"], kw["gausssigma"],
			extract<std::string>(kw["wind"]), kw["windparams"]);
	const auto irr = make_neutron_star_self_irradiation_arguments(
			extract<double>(kw["Cirr"]),
			extract<double>(kw["irrindex"]),
			extract<double>(kw["Cirrcold"]),
			extract<double>(kw["irrindexcold"]),
			extract<double>(kw["h2rcold"]),
			extract<std::string>(kw["angulardistdisk"]),
			extract<std::string>(kw["angulardistns"]));
	const auto flux = make_flux_arguments(
			extract<double>(kw["colourfactor"]),
			extract<double>(kw["emin"]), extract<double>(kw["emax"]),
			extract<double>(kw["staralbedo"]),
			extract<double>(kw["inclination"]), extract<double>(kw["ephemerist0"]), extract<double>(kw["distance"]));
	const auto calc = make_calculation_arguments(
			extract<double>(kw["inittime"]),
			extract<double>(kw["time"]), kw["tau"],
			extract<unsigned int>(kw["Nx"]), extract<std::string>(kw["gridscale"]),
			extract<unsigned short>(kw["starlod"]),
			kw["eps"]);
	return make_freddi_neutron_star_arguments(*general, *basic, *disk, *irr, *flux, *calc, *ns_args);
}

object raw_make_neutron_star_evolution(const tuple& args, const dict& kwargs) {
	check_args(args);

	// It should be static, but segfault
	const auto kwdefaults = neutron_star_evolution_kwdefaults();
	check_kwargs(args, kwargs, neutron_star_evolution_required_args(), kwdefaults);

	dict kw;
	kw.update(kwdefaults);
	kw.update(kwargs);

	const auto freddi_args = make_freddi_neutron_star_arguments(kw);

	args[0].attr("_kwargs") = kw;
	return args[0].attr("__init__")(*freddi_args);
}


double (FreddiNeutronStarEvolution::*fp_getter)() const = &FreddiNeutronStarEvolution::fp;
double (FreddiNeutronStarEvolution::*eta_ns_getter)() const = &FreddiNeutronStarEvolution::eta_ns;


void wrap_evolution() {
	class_<FreddiEvolution, bases<FreddiState> >("_Freddi", no_init)
		.def("__init__", raw_function(&raw_make_evolution))
		.def(init<const FreddiArguments&>())
		.def("__iter__", iterator<FreddiEvolution>())
		.def("_required_args", evolution_required_args, "Mock values for non-scientific calls")
		.staticmethod("_required_args")
	;
	class_<FreddiNeutronStarEvolution, bases<FreddiEvolution> >("_FreddiNeutronStar", no_init)
		.def("__init__", raw_function(&raw_make_neutron_star_evolution))
		.def(init<const FreddiNeutronStarArguments&>())
		.def("__iter__", iterator<FreddiNeutronStarEvolution>())
		.def("_required_args", neutron_star_evolution_required_args, "Mock values for non-scientific calls")
		.staticmethod("_required_args")
		.add_property("mu_magn", &FreddiNeutronStarEvolution::mu_magn)
		.add_property("R_cor", &FreddiNeutronStarEvolution::R_cor)
		.add_property("R_dead", &FreddiNeutronStarEvolution::R_dead)
		.add_property("inverse_beta", &FreddiNeutronStarEvolution::inverse_beta)
		.add_property("Fmagn", make_function(&FreddiNeutronStarEvolution::Fmagn, return_value_policy<copy_const_reference>()))
		.add_property("dFmagn_dh", make_function(&FreddiNeutronStarEvolution::dFmagn_dh, return_value_policy<copy_const_reference>()))
		.add_property("d2Fmagn_dh2", make_function(&FreddiNeutronStarEvolution::d2Fmagn_dh2, return_value_policy<copy_const_reference>()))
		.add_property("eta_ns", eta_ns_getter)
		.add_property("fp", fp_getter)
		.add_property("T_hot_spot", &FreddiNeutronStarEvolution::T_hot_spot)
		.add_property("Lbol_ns", &FreddiNeutronStarEvolution::Lbol_ns)
		.add_property("Lx_ns", &FreddiNeutronStarEvolution::Lx_ns)
	;
}
