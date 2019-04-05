#include <iostream> // TODO: remove
#include <sstream>
#include <string>

#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>

#include <freddi_evolution.hpp>
#include <unit_transformation.hpp>

#include "arguments.hpp"
#include "freddi_evolution.hpp"

using namespace boost::python;


dict evolution_kwdefaults() {
	dict kw;

	kw["__cgs"] = true;

	kw["alpha"] = BasicDiskBinaryArguments::default_alpha;
	kw["Mx"] = BasicDiskBinaryArguments::default_Mx;
	kw["kerr"] = BasicDiskBinaryArguments::default_kerr;
	kw["Mopt"] = BasicDiskBinaryArguments::default_Mopt;
	kw["period"] = BasicDiskBinaryArguments::default_period;
	kw["rin"] = object();
	kw["rout"] = object();

	kw["opacity"] = DiskStructureArguments::default_opacity;
	kw["Mdotout"] = DiskStructureArguments::default_Mdotout;
	kw["boundcond"] = DiskStructureArguments::default_boundcond;
	kw["Thot"] = DiskStructureArguments::default_Thot;
	kw["initialcond"] = DiskStructureArguments::default_initialcond;
	kw["F0"] = DiskStructureArguments::default_F0;
	kw["powerorder"] = DiskStructureArguments::default_powerorder;
	kw["gaussmu"] = DiskStructureArguments::default_gaussmu;
	kw["gausssigma"] = DiskStructureArguments::default_gausssigma;
	kw["Mdisk0"] = object();
	kw["Mdot0"] = object();
	kw["wind"] = DiskStructureArguments::default_wind;
	kw["windparams"] = tuple();

	kw["Cirr"] = SelfIrradiationArguments::default_Cirr;
	kw["irrfactortype"] = SelfIrradiationArguments::default_irrfactortype;

	kw["colourfactor"] = FluxArguments::default_colourfactor;
	kw["emin"] = FluxArguments::default_emin;
	kw["emax"] = FluxArguments::default_emax;
	kw["inclination"] = FluxArguments::default_inclination;
	kw["distance"] = FluxArguments::default_distance;
	kw["lambdas"] = tuple();

	kw["time"] = CalculationArguments::default_time;
	kw["tau"] = CalculationArguments::default_tau;
	kw["Nx"] = CalculationArguments::default_Nx;
	kw["gridscale"] = CalculationArguments::default_gridscale;
	kw["eps"] = object();

	return kw;
}


void check_args(const tuple &args) {
	if (len(args) > 1) {
		PyErr_SetString(PyExc_TypeError, "Positional arguments are denied");
		throw error_already_set();
	}
}


void check_kwargs(const tuple& args, const dict& kwargs, const dict& kwdefaults) {
	stl_input_iterator<object> begin(kwargs), end;
	for (auto key = begin; key != end; ++key) {
		if (! kwdefaults.has_key(*key)) {
			std::stringstream msg;
			msg << extract<char*>(str(args[0])) << " got an unexpected keyword argument " << extract<char*>(*key);
			PyErr_SetString(PyExc_TypeError, msg.str().data());
			throw error_already_set();
		}
	}
}


boost::shared_ptr<FreddiArguments> make_freddi_arguments(dict& kw) {
	if (! extract<bool>(kw["__cgs"])) {
		kw["Mx"] = sunToGram(extract<double>(kw["Mx"]));
		kw["Mopt"] = sunToGram(extract<double>(kw["Mopt"]));
		kw["period"] = dayToS(extract<double>(kw["period"]));
		if (object(kw["rin"]).ptr() != object().ptr()) {
			kw["rin"] = rgToCm(extract<double>(kw["rin"]), extract<double>(kw["Mx"]));
		}
		if (object(kw["rout"]).ptr() != object().ptr()) {
			kw["rout"] = sunToCm(extract<double>(kw["rout"]));
		}
		kw["emin"] = kevToHertz(extract<double>(kw["emin"]));
		kw["emax"] = kevToHertz(extract<double>(kw["emax"]));
		kw["distance"] = kpcToCm(extract<double>(kw["distance"]));
		for (size_t i = 0; i < len(kw["lambdas"]); ++i) {
			kw["lambdas"][i] = angstromToCm(extract<double>(kw["lambdas"][i]));
		}
		kw["time"] = dayToS(extract<double>(kw["time"]));
		kw["tau"] = dayToS(extract<double>(kw["tau"]));
	}

	const auto general = make_general_arguments();
	const auto basic = make_basic_disk_binary_arguments(
			extract<double>(kw["alpha"]),
			extract<double>(kw["Mx"]), extract<double>(kw["kerr"]),
			extract<double>(kw["Mopt"]), extract<double>(kw["period"]),
			kw["rin"], kw["rout"]);
	const auto disk = make_disk_structure_arguments(
			*basic,
			extract<std::string>(kw["opacity"]),
			extract<double>(kw["Mdotout"]),
			extract<std::string>(kw["boundcond"]), extract<double>(kw["Thot"]),
			extract<std::string>(kw["initialcond"]), extract<double>(kw["F0"]),
			extract<double>(kw["powerorder"]), extract<double>(kw["gaussmu"]), extract<double>(kw["gausssigma"]),
			kw["Mdisk0"], kw["Mdot0"],
			extract<std::string>(kw["wind"]), kw["windparams"]);
	const auto irr = make_self_irradiation_arguments(
			extract<double>(kw["Cirr"]),
			extract<std::string>(kw["irrfactortype"]));
	const auto flux = make_flux_arguments(
			extract<double>(kw["colourfactor"]),
			extract<double>(kw["emin"]), extract<double>(kw["emax"]),
			extract<double>(kw["inclination"]), extract<double>(kw["distance"]),
			kw["lambdas"]);
	const auto calc = make_calculation_arguments(
			extract<double>(kw["time"]), extract<double>(kw["tau"]),
			extract<unsigned int>(kw["Nx"]), extract<std::string>(kw["gridscale"]),
			kw["eps"]);
	// To avoid copy, create FreddiArguments constructor that accepts shared_ptr
	return make_freddi_arguments(*general, *basic, *disk, *irr, *flux, *calc);
}


object raw_make_evolution(const tuple& args, const dict& kwargs) {
	check_args(args);

	// It should be static, but segfault
	const auto kwdefaults = evolution_kwdefaults();
	check_kwargs(args, kwargs, kwdefaults);

	dict kw;
	kw.update(kwdefaults);
	kw.update(kwargs);

	const auto freddi_args = make_freddi_arguments(kw);

	args[0].attr("_kwargs") = kw;
	return args[0].attr("__init__")(*freddi_args);
}


dict neutron_star_evolution_kwdefaults() {
	dict kw;
	kw.update(evolution_kwdefaults());

	kw["Rx"] = NeutronStarArguments::default_Rx;
	kw["freqx"] = NeutronStarArguments::default_freqx;
	kw["Bx"] = NeutronStarArguments::default_Bx;
	kw["hotspotarea"] = NeutronStarArguments::default_hotspotarea;
	kw["epsilonAlfven"] = NeutronStarArguments::default_epsilonAlfven;
	kw["inversebeta"] = NeutronStarArguments::default_inversebeta;
	kw["Rdead"] = NeutronStarArguments::default_Rdead;
	kw["fptype"] = NeutronStarArguments::default_fptype;

	return kw;
}

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(dict& kw) {
	const auto freddi_args = make_freddi_arguments(kw);
	const auto ns_args = make_neutron_star_arguments(
			extract<double>(kw["Rx"]), extract<double>(kw["freqx"]), extract<double>(kw["Bx"]), extract<double>(kw["hotspotarea"]),
			extract<double>(kw["epsilonAlfven"]), extract<double>(kw["inversebeta"]), extract<double>(kw["Rdead"]),
			extract<std::string>(kw["fptype"]));
	return make_freddi_neutron_star_arguments(*freddi_args, *ns_args);
}

object raw_make_neutron_star_evolution(const tuple& args, const dict& kwargs) {
	check_args(args);

	// It should be static, but segfault
	const auto kwdefaults = neutron_star_evolution_kwdefaults();
	check_kwargs(args, kwargs, kwdefaults);

	dict kw;
	kw.update(kwdefaults);
	kw.update(kwargs);

	const auto freddi_args = make_freddi_neutron_star_arguments(kw);

	args[0].attr("_kwargs") = kw;
	return args[0].attr("__init__")(*freddi_args);
}


double (FreddiNeutronStarEvolution::*fp_getter)() const = &FreddiNeutronStarEvolution::fp;


void wrap_evolution() {
	class_<FreddiEvolution, bases<FreddiState> >("_Freddi", no_init)
		.def("__init__", raw_function(&raw_make_evolution))
		.def(init<const FreddiArguments&>())
		.def("__iter__", iterator<FreddiEvolution>())
	;
	class_<FreddiNeutronStarEvolution, bases<FreddiEvolution> >("_FreddiNeutronStar", no_init)
		.def("__init__", raw_function(&raw_make_neutron_star_evolution))
		.def(init<const FreddiNeutronStarArguments&>())
		.def("__iter__", iterator<FreddiNeutronStarEvolution>())
		.add_property("mu_magn", &FreddiNeutronStarEvolution::mu_magn)
		.add_property("R_dead", &FreddiNeutronStarEvolution::R_dead)
		.add_property("F_dead", &FreddiNeutronStarEvolution::F_dead)
		.add_property("R_cor", &FreddiNeutronStarEvolution::R_cor)
		.add_property("inverse_beta", &FreddiNeutronStarEvolution::inverse_beta)
		.add_property("Fmagn", make_function(&FreddiNeutronStarEvolution::Fmagn, return_value_policy<copy_const_reference>()))
		.add_property("dFmagn_dh", make_function(&FreddiNeutronStarEvolution::dFmagn_dh, return_value_policy<copy_const_reference>()))
		.add_property("d2Fmagn_dh2", make_function(&FreddiNeutronStarEvolution::d2Fmagn_dh2, return_value_policy<copy_const_reference>()))
		.add_property("eta_ns", &FreddiNeutronStarEvolution::eta_ns)
		.add_property("fp", fp_getter)
		.add_property("T_hot_spot", &FreddiNeutronStarEvolution::T_hot_spot)
		.add_property("Lbol_ns", &FreddiNeutronStarEvolution::Lbol_ns)
		.add_property("Lx_ns", &FreddiNeutronStarEvolution::Lx_ns)
	;
}
