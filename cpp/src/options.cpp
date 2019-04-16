#include <algorithm>  // transform

#include "options.hpp"

namespace po = boost::program_options;


GeneralOptions::GeneralOptions(const po::variables_map& vm):
		GeneralArguments(
				vm["prefix"].as<std::string>(),
				vm["dir"].as<std::string>(),
				(vm.count("fulldata") > 0)) {}

po::options_description GeneralOptions::description() {
	po::options_description od("General options");
	od.add_options()
			( "help,h", "Produce help message" )
			( "prefix", po::value<std::string>()->default_value(default_prefix), "Set prefix for output filenames. Output file with distribution of parameters over time is PREFIX.dat" )
			( "dir,d", po::value<std::string>()->default_value(default_dir), "Choose the directory to write output files. It should exist" )
			( "fulldata", "Output files PREFIX_%d.dat with radial structure for every time step. Default is to output only PREFIX.dat with global disk parameters for every time step" )
			;
	return od;
}


BasicDiskBinaryOptions::BasicDiskBinaryOptions(const po::variables_map &vm):
		BasicDiskBinaryArguments(
				vm["alpha"].as<double>(),
				sunToGram(vm["Mx"].as<double>()),
				vm["kerr"].as<double>(),
				sunToGram(vm["Mopt"].as<double>()),
				dayToS(vm["period"].as<double>()),
				rinInitializer(vm),
				routInitializer(vm)) {
	if (rin >= rout) {
		throw po::invalid_option_value("rin should be smaller rout");
	}
}

double BasicDiskBinaryOptions::rinInitializer(const po::variables_map &vm) {
	const double Mx = sunToGram(vm["Mx"].as<double>());
	if (vm.count("rin")) {
		return rgToCm(vm["rin"].as<double>(), Mx);
	}
	const double kerr = vm["kerr"].as<double>();
	return rinFromMxKerr(Mx, kerr);
}

double BasicDiskBinaryOptions::routInitializer(const po::variables_map &vm) {
	if (vm.count("rout")) {
		return sunToCm(vm["rout"].as<double>());
	}
	const double Mx = sunToGram(vm["Mx"].as<double>());
	const double Mopt = sunToGram(vm["Mopt"].as<double>());
	const double period = dayToS(vm["period"].as<double>());
	return routFromMxMoptPeriod(Mx, Mopt, period);
}

po::options_description BasicDiskBinaryOptions::description() {
	po::options_description od("Basic binary and disk parameter");
	od.add_options()
			( "alpha,a", po::value<double>()->default_value(default_alpha), "Alpha parameter of Shakura-Sunyaev model" )
			( "Mx,M", po::value<double>()->default_value(gramToSun(default_Mx)), "Mass of the central object, in the units of solar masses" )
			( "kerr", po::value<double>()->default_value(default_kerr), "Dimensionless Kerr parameter of the black hole" )
			( "Mopt",	po::value<double>()->default_value(gramToSun(default_Mopt)), "Mass of the optical star, in units of solar masses" )
			( "period,P", po::value<double>()->default_value(sToDay(default_period)), "Orbital period of the binary system, in units of days" )
			( "rin", po::value<double>(), "Inner radius of the disk, in the units of the Schwarzschild radius of the central object 2GM/c^2. If it isn't set then the radius of ISCO orbit is used defined by --Mx and --kerr values" )
			( "rout,R", po::value<double>(), "Outer radius of the disk, in units of solar radius. If it isn't set then the tidal radius is used defined by --Mx, --Mopt and --period values" )
			;
	return od;
}


DiskStructureOptions::DiskStructureOptions(const po::variables_map &vm, const BasicDiskBinaryArguments& bdb_args):
		DiskStructureArguments(
				bdb_args,
				vm["opacity"].as<std::string>(),
				vm["Mdotout"].as<double>(),
				vm["boundcond"].as<std::string>(),
				vm["Thot"].as<double>(),
				vm["initialcond"].as<std::string>(),
				vm["F0"].as<double>(),
				vm["powerorder"].as<double>(),
				vm["gaussmu"].as<double>(),
				vm["gausssigma"].as<double>(),
				(vm.count("Mdisk0") > 0),
				(vm.count("Mdot0") > 0),
				Mdisk0Initializer(vm),
				Mdot0Initializer(vm),
				"no", {}) {}  // wind

double DiskStructureOptions::Mdisk0Initializer(const po::variables_map& vm) {
	if (vm.count("Mdisk0") > 0) {
		return vm["Mdisk0"].as<double>();
	}
	return -1;
}

double DiskStructureOptions::Mdot0Initializer(const po::variables_map& vm) {
	if (vm.count("Mdot0") > 0) {
		return vm["Mdot0"].as<double>();
	}
	return -1;
}

po::options_description DiskStructureOptions::description() {
	po::options_description od("Parameters of the disk mode");
	od.add_options()
			( "opacity,O", po::value<std::string>()->default_value(default_opacity), "Opacity law: Kramers (varkappa ~ rho / T^7/2) or OPAL (varkappa ~ rho / T^5/2)" )
			( "Mdotout", po::value<double>()->default_value(default_Mdotout), "Accretion rate onto the disk through its outer radius" )
			( "boundcond", po::value<std::string>()->default_value(default_boundcond), "Outer boundary movement condition\n\n"
																					   "Values:\n"
																					   "  Teff: outer radius of the disk moves inwards to keep photosphere temperature of the disk larger than some value. This value is specified by --Thot option\n"
																					   "  Tirr: outer radius of the disk moves inwards to keep irradiation flux of the disk larger than some value. The value of this minimal irradiation flux is [Stefan-Boltzmann constant] * Tirr^4, where Tirr is specified by --Thot option" ) // fourSigmaCrit, MdotOut
			( "Thot", po::value<double>()->default_value(default_Thot), "Minimum photosphere or irradiation temperature at the outer edge of the hot disk, Kelvin. For details see --boundcond description" )
			( "initialcond", po::value<std::string>()->default_value(default_initialcond), "Type of the initial condition for viscous torque F or surface density Sigma\n\n"
																						   "Values:\n"
																						   "  powerF: F ~ xi^powerorder, powerorder is specified by --powerorder option\n" // power option does the same
																						   "  powerSigma: Sigma ~ xi^powerorder, powerorder is specified by --powerorder option\n"
																						   "  sinusF: F ~ sin( xi * pi/2 )\n" // sinus option does the same
																						   "  gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2), mu and sigma are specified by --gaussmu and --gausssigma options\n"
																						   "  quasistat: F ~ f(h/h_out) * xi * h_out/h, where f is quasi-stationary solution found in Lipunova & Shakura 2000. f(xi=0) = 0, df/dxi(xi=1) = 0\n\n"
																						   "Here xi is (h - h_in) / (h_out - h_in)\n") // sinusparabola, sinusgauss
			( "F0", po::value<double>()->default_value(default_F0), "Initial maximum viscous torque in the disk, dyn*cm. Can be overwritten via --Mdisk0 and --Mdot0" )
			( "Mdisk0", po::value<double>(), "Initial disk mass, g. If both --F0 and --Mdisk0 are specified then --Mdisk0 is used. If both --Mdot0 and --Mdisk0 are specified then --Mdot0 is used" )
			( "Mdot0", po::value<double>(), "Initial mass accretion rate through the inner radius, g/s. If --F0, --Mdisk0 and --Mdot0 are specified then --Mdot0 is used. Works only when --initialcond is set to sinusF or quasistat" )
			( "powerorder", po::value<double>()->default_value(default_powerorder), "Parameter for the powerlaw initial condition distribution. This option works only with --initialcond=powerF or powerSigma" )
			( "gaussmu", po::value<double>()->default_value(default_gaussmu), "Position of the maximum for Gauss distribution, positive number not greater than unity. This option works only with --initialcond=gaussF" )
			( "gausssigma", po::value<double>()->default_value(default_gausssigma), "Width of for Gauss distribution. This option works only with --initialcond=gaussF" )
			;
	return od;
}


SelfIrradiationOptions::SelfIrradiationOptions(const po::variables_map &vm, const DiskStructureArguments &dsa_args):
		SelfIrradiationArguments(
				vm["Cirr"].as<double>(),
				vm["irrfactortype"].as<std::string>()) {
	if (Cirr <= 0. && dsa_args.boundcond == "Tirr") {
		throw po::error("It is obvious to use nonpositive --Cirr with --boundcond=Tirr");
	}
	if (irrfactortype != "const" && irrfactortype != "square") {
		throw po::invalid_option_value("--irrfactortype has invalid value");
	}
}

po::options_description SelfIrradiationOptions::description() {
	po::options_description od("Parameters of self-irradiatio");
	od.add_options()
			( "Cirr", po::value<double>()->default_value(default_Cirr), "Irradiation factor" )
			( "irrfactortype", po::value<std::string>()->default_value(default_irrfactortype), "Type of irradiation factor Cirr\n\n"
																							   "Values:\n"
																							   "  const: doesn't depend on disk shape:\n[rad. flux] = Cirr  L / (4 pi r^2)\n"
																							   "  square: Cirr depends on the disk relative half-thickness:\n[rad. flux] = Cirr (z/r)^2 L / (4 pi r^2)\n\n"
																							   "Here L is bolometric Luminosity:\nL = eta Mdot c^2" )
			;
	return od;
}


FluxOptions::FluxOptions(const po::variables_map &vm):
		FluxArguments(
				vm["colourfactor"].as<double>(),
				kevToHertz(vm["emin"].as<double>()),
				kevToHertz(vm["emax"].as<double>()),
				vm["inclination"].as<double>(),
				kpcToCm(vm["distance"].as<double>()),
				lambdasInitializer(vm)) {}

vecd FluxOptions::lambdasInitializer(const po::variables_map &vm) const {
	if (vm.count("lambda") > 0) {
		vecd lambdas(vm["lambda"].as<vecd>());
		transform(lambdas.begin(), lambdas.end(), lambdas.begin(), angstromToCm);
		return lambdas;
	}
	return vecd();
}

po::options_description FluxOptions::description() {
	po::options_description od("Parameters of flux calculation");
	od.add_options()
			( "colourfactor", po::value<double>()->default_value(default_colourfactor), "Colour factor to calculate X-ray flux"  )
			( "emin", po::value<double>()->default_value(hertzToKev(default_emin)), "Minimum energy of X-ray band, keV" )
			( "emax", po::value<double>()->default_value(hertzToKev(default_emax)), "Maximum energy of X-ray band, keV" )
			( "inclination,i", po::value<double>()->default_value(default_inclination), "Inclination of the system, degrees" )
			( "distance", po::value<double>()->default_value(cmToKpc(default_distance)), "Distance to the system, kpc" )
			( "lambda", po::value<vecd>()->multitoken(), "Wavelength to calculate Fnu, Angstrom. You can use this option multiple times. For each lambda one additional column with values of spectral flux density Fnu [erg/s/cm^2/Hz] is produced" )
			;
	return od;
}


CalculationOptions::CalculationOptions(const po::variables_map &vm):
		CalculationArguments(
				dayToS(vm["time"].as<double>()),
				dayToS(vm["tau"].as<double>()),
				vm["Nx"].as<unsigned int>(),
				vm["gridscale"].as<std::string>()) {
	if (gridscale != "log" && gridscale != "linear") {
		throw po::invalid_option_value("Invalid --gridscale value");
	}
}

po::options_description CalculationOptions::description() {
	po::options_description od("Parameters of disk evolution calculation");
	od.add_options()
			( "time,T", po::value<double>()->default_value(sToDay(default_time)), "Time interval to calculate evolution, days" )
			( "tau",	po::value<double>()->default_value(sToDay(default_tau)), "Time step, days" )
			( "Nx",	po::value<unsigned int>()->default_value(default_Nx), "Size of calculation grid" )
			( "gridscale", po::value<std::string>()->default_value(default_gridscale), "Type of grid for angular momentum h: log or linear" )
			;
	return od;
}



FreddiOptions::FreddiOptions(const po::variables_map& vm) {
	general.reset(new GeneralOptions(vm));
	basic.reset(new BasicDiskBinaryOptions(vm));
	disk.reset(new DiskStructureOptions(vm, *basic));
	irr.reset(new SelfIrradiationOptions(vm, *disk));
	flux.reset(new FluxOptions(vm));
	calc.reset(new CalculationOptions(vm));
}

po::options_description FreddiOptions::description() {
	po::options_description desc("Freddi: numerical calculation of accretion disk evolution");
	desc.add(GeneralOptions::description());
	desc.add(BasicDiskBinaryOptions::description());
	desc.add(DiskStructureOptions::description());
	desc.add(SelfIrradiationOptions::description());
	desc.add(FluxOptions::description());
	desc.add(CalculationOptions::description());
	return desc;
}


NeutronStarOptions::NeutronStarOptions(const po::variables_map &vm):
		NeutronStarArguments(
				vm["Rx"].as<double>(),
				vm["freqx"].as<double>(),
				vm["Bx"].as<double>(),
				vm["hotspotarea"].as<double>(),
				vm["epsilonAlfven"].as<double>(),
				vm["inversebeta"].as<double>(),
				vm["Rdead"].as<double>(),
				vm["fptype"].as<std::string>(),
				{}) {}

po::options_description NeutronStarOptions::description() {
	po::options_description od("Parameters of accreting neutron star");
	od.add_options()
			( "Rx", po::value<double>()->default_value(default_Rx), "Accretor radius, cm" )
			( "freqx", po::value<double>()->default_value(default_freqx), "Accretor rotation frequency, Hz. This parameter is not linked to --kerr, agree them yourself" )
			( "Bx", po::value<double>()->default_value(default_Bx), "Accretor polar magnetic induction, G" )
			( "hotspotarea", po::value<double>()->default_value(default_hotspotarea), "??? Relative area of hot spot on the accretor" )
			( "epsilonAlfven", po::value<double>()->default_value(default_epsilonAlfven), "Factor in Alfven radius formula" )
			( "inversebeta", po::value<double>()->default_value(default_inversebeta), "???" )
			( "Rdead", po::value<double>()->default_value(default_Rdead), "Maximum inner radius of the disk that can be obtained, it characterises minimum torque in the dead disk, cm" )
			( "fptype", po::value<std::string>()->default_value(default_fptype), "??? Accretor Mdot fraction" )
			;
	return od;
}


FreddiNeutronStarOptions::FreddiNeutronStarOptions(const po::variables_map &vm) {
	general.reset(new GeneralOptions(vm));
	basic.reset(new BasicDiskBinaryOptions(vm));
	disk.reset(new DiskStructureOptions(vm, *basic));
	irr.reset(new SelfIrradiationOptions(vm, *disk));
	flux.reset(new FluxOptions(vm));
	calc.reset(new CalculationOptions(vm));
	ns.reset(new NeutronStarOptions(vm));
}

po::options_description FreddiNeutronStarOptions::description() {
	po::options_description desc("Freddi NS: numerical calculation of accretion disk evolution");
	desc.add(GeneralOptions::description());
	desc.add(BasicDiskBinaryOptions::description());
	desc.add(DiskStructureOptions::description());
	desc.add(NeutronStarOptions::description());
	desc.add(SelfIrradiationOptions::description());
	desc.add(FluxOptions::description());
	desc.add(CalculationOptions::description());
	return desc;
}

