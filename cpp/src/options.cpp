#include <algorithm>  // transform
#include <vector>

#include <boost/algorithm/string.hpp> // split is_any_of

#include "options.hpp"

namespace po = boost::program_options;


GeneralOptions::GeneralOptions(const po::variables_map& vm):
		GeneralArguments(
				vm["prefix"].as<std::string>(),
				vm["dir"].as<std::string>(),
				vm["precision"].as<unsigned int>(),
				vm["tempsparsity"].as<unsigned int>(),
				(vm.count("fulldata") > 0),
				(vm.count("stdout") > 0)) {}

po::options_description GeneralOptions::description() {
	po::options_description od("General options");
	od.add_options()
			( "help,h", "Produce help message" )
			( "config", po::value<std::string>(), "Set additional configuration filepath" )
			( "prefix", po::value<std::string>()->default_value(default_prefix), "Set prefix for output filenames. Output file with distribution of parameters over time is PREFIX.dat" )
			( "stdout", "Output temporal distribution to stdout instead of PREFIX.dat file" )
			( "dir,d", po::value<std::string>()->default_value(default_dir), "Choose the directory to write output files. It should exist" )
			( "precision", po::value<unsigned int>()->default_value(default_output_precision), "Number of digits to print into output files" )
			( "tempsparsity", po::value<unsigned int>()->default_value(default_temp_sparsity_output), "Output every k-th time moment" )
			( "fulldata", "Output files PREFIX_%d.dat with radial structure for every time step. Default is to output only PREFIX.dat with global disk parameters for every time step" )
			;
	return od;
}


BasicDiskBinaryOptions::BasicDiskBinaryOptions(const po::variables_map &vm):
		BasicDiskBinaryArguments(
				vm["alpha"].as<double>(),
				varToOpt<double>(vm, "alphacold"),
				sunToGram(vm["Mx"].as<double>()),
				vm["kerr"].as<double>(),
				dayToS(vm["period"].as<double>()),
				sunToGram(vm["Mopt"].as<double>()),
				vm["rochelobefill"].as<double>(),
				vm["Topt"].as<double>(),
				rinInitializer(vm),
				routInitializer(vm),
				riscoInitializer(vm)) {
	if (rin >= rout) {
		throw po::invalid_option_value("rin should be smaller rout");
	}
}

std::optional<double> BasicDiskBinaryOptions::rinInitializer(const po::variables_map &vm) {
	if (vm.count("rin")) {
		const double Mx = sunToGram(vm["Mx"].as<double>());
		return rgToCm(vm["rin"].as<double>(), Mx);
	}
	return {};
}

std::optional<double> BasicDiskBinaryOptions::routInitializer(const po::variables_map &vm) {
	if (vm.count("rout")) {
		return sunToCm(vm["rout"].as<double>());
	}
	return {};
}

std::optional<double> BasicDiskBinaryOptions::riscoInitializer(const po::variables_map &vm) {
	if (vm.count("risco")) {
		const double Mx = sunToGram(vm["Mx"].as<double>());
		return rgToCm(vm["risco"].as<double>(), Mx);
	}
	return {};
}

po::options_description BasicDiskBinaryOptions::description() {
	po::options_description od("Basic binary and disk parameter");
	od.add_options()
			( "alpha,a", po::value<double>()->required(), "Alpha parameter of Shakura-Sunyaev model" )
			( "alphacold", po::value<double>(), "Alpha parameter of cold disk, currently it is used only for Sigma_minus, see --Qirr2Qvishot. Default is --alpha values divided by ten" )  // default_alpha_to_alphacold
			( "Mx,M", po::value<double>()->required(), "Mass of the central object, in the units of solar masses" )
			( "kerr", po::value<double>()->default_value(default_kerr), "Dimensionless Kerr parameter of the black hole" )
			( "Mopt",	po::value<double>()->required(), "Mass of the optical star, in units of solar masses" )
			( "rochelobefill", po::value<double>()->default_value(default_roche_lobe_fill), "Dimensionless factor describing a size of the optical star. Polar radius of the star is rochelobefill * (polar radius of critical Roche lobe)" )
			( "Topt", po::value<double>()->default_value(default_Topt), "Thermal temperature of the optical star, in units of kelvins" )
			( "period,P", po::value<double>()->required(), "Orbital period of the binary system, in units of days" )
			( "rin", po::value<double>(), "Inner radius of the disk, in the units of the gravitational radius of the central object GM/c^2. If it isn't set then the radius of ISCO orbit is used defined by --Mx and --kerr values" )
			( "rout,R", po::value<double>(), "Outer radius of the disk, in units of solar radius. If it isn't set then the tidal radius is used defined by --Mx, --Mopt and --period values" )
			( "risco", po::value<double>(), "Innermost stable circular orbit, in units of gravitational radius of the central object GM/c^2. If it isn't set then the radius of ISCO orbit is used defined by --Mx and --kerr values" )
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
				std::pow(vm["Qirr2Qvishot"].as<double>(), 0.25),
				vm["initialcond"].as<std::string>(),
				varToOpt<double>(vm, "F0"),
				varToOpt<double>(vm, "Mdisk0"),
				varToOpt<double>(vm, "Mdot0"),
				varToOpt<double>(vm, "powerorder"),
				varToOpt<double>(vm, "gaussmu"),
				varToOpt<double>(vm, "gausssigma"),
				vm["windtype"].as<std::string>(),
				windparamsInitializer(vm)) {}


pard DiskStructureOptions::windparamsInitializer(const po::variables_map& vm) {
	const auto windtype = vm["windtype"].as<std::string>();

	if (windtype == "no") {
		return {};
	}
	if (windtype == "SS73C") {
		return {};
	}
	if (windtype == "Woods1996"){
		if (vm.count("windXi") == 0) {
			throw po::error("--windXi is required if --windtype=Woods1996");
		}
		if (vm.count("windTic") == 0) {
			throw po::error("--windTic is required if --windtype=Woods1996");
		}
		if (vm.count("windPow") == 0) {
			throw po::error("--windPow is required if --windtype=Woods1996");
		}
		return {
				{"Xi_max", vm["windXi"].as<double>()},
				{"T_iC", vm["windTic"].as<double>()},
				{"W_pow", vm["windPow"].as<double>()}
		};
	}
	if (windtype == "toy"){
		if (vm.count("windPow") == 0) {
			throw po::error("--windPow is required if --windtype=toy");
		}
		return {

				{"W_pow", vm["windPow"].as<double>()}
		};
	}
	

	throw po::invalid_option_value("Unknown --windtype=" + windtype);
}


po::options_description DiskStructureOptions::description() {
	po::options_description od("Parameters of the disk mode");
	od.add_options()
			( "opacity,O", po::value<std::string>()->default_value(default_opacity), "Opacity law: Kramers (varkappa ~ rho / T^7/2) or OPAL (varkappa ~ rho / T^5/2)" )
			( "Mdotout", po::value<double>()->default_value(default_Mdotout), "Accretion rate onto the disk through its outer radius" )
			( "boundcond", po::value<std::string>()->default_value(default_boundcond),
					"Outer boundary movement condition\n\n"
					"Values:\n"
					"  Teff: outer radius of the disk moves inwards to keep photosphere temperature of the disk larger than some value. This value is specified by --Thot option\n"
					"  Tirr: outer radius of the disk moves inwards to keep irradiation flux of the disk larger than some value. The value of this minimal irradiation flux is [Stefan-Boltzmann constant] * Tirr^4, where Tirr is specified by --Thot option" ) // fourSigmaCrit, MdotOut
			( "Thot", po::value<double>()->default_value(default_Thot), "Minimum photosphere or irradiation temperature at the outer edge of the hot disk, Kelvin. For details see --boundcond description" )
			( "Qirr2Qvishot", po::value<double>()->default_value(m::pow<4>(default_Tirr2Tvishot)), "Minimum Qirr / Qvis ratio at the outer edge of the hot disk to switch evolution from temperature-based regime to Sigma_minus-based regime (see Eq. A.1 in Lasota et al. 2008, --alphacold value is used as alpha parameter)" )
			( "initialcond", po::value<std::string>()->default_value(default_initialcond),
					"Type of the initial condition for viscous torque F or surface density Sigma\n\n"
					"Values:\n"
					"  powerF: F ~ xi^powerorder, powerorder is specified by --powerorder option\n" // power does the same
					"  linearF: F ~ xi, specific case of powerF but can be normalised by --Mdot0, see its description for details" // linear does the same
					"  powerSigma: Sigma ~ xi^powerorder, powerorder is specified by --powerorder option\n"
					"  sineF: F ~ sin( xi * pi/2 )\n" // sinus option does the same
					"  gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2), mu and sigma are specified by --gaussmu and --gausssigma options\n"
					"  quasistat: F ~ f(h/h_out) * xi * h_out/h, where f is quasi-stationary solution found in Lipunova & Shakura 2000. f(xi=0) = 0, df/dxi(xi=1) = 0\n\n"
					"Here xi is (h - h_in) / (h_out - h_in)\n")
			( "F0", po::value<double>(), "Initial maximum viscous torque in the disk, dyn*cm. Can be overwritten via --Mdisk0 and --Mdot0" )
			( "Mdisk0", po::value<double>(), "Initial disk mass, g. If both --F0 and --Mdisk0 are specified then --Mdisk0 is used. If both --Mdot0 and --Mdisk0 are specified then --Mdot0 is used" )
			( "Mdot0", po::value<double>(), "Initial mass accretion rate through the inner radius, g/s. If --F0, --Mdisk0 and --Mdot0 are specified then --Mdot0 is used. Works only when --initialcond is set to linearF, sinusF or quasistat" )
			( "powerorder", po::value<double>(), "Parameter for the powerlaw initial condition distribution. This option works only with --initialcond=powerF or powerSigma" )
			( "gaussmu", po::value<double>(), "Position of the maximum for Gauss distribution, positive number not greater than unity. This option works only with --initialcond=gaussF" )
			( "gausssigma", po::value<double>(), "Width of for Gauss distribution. This option works only with --initialcond=gaussF" )
			( "windtype", po::value<std::string>()->default_value(default_wind),
			        "Type of the wind\n\n"
					"Values:\n"
					"  no: no wind\n"
					"  SS73C: super-Eddington spherical wind from Shakura-Sunyaev 1973\n"
					"  Janiuk2015: super-Eddington Janiuk et al. 2015\n"
					"  Woods1996: thermal wind Woods et al. 1996. Requires --windXi, --windTic and --windPow to be specified"
					"  toy: a toy wind model used in arXiv:2105.11974, the mass loss rate is proportional to the central accretion rate. Requires --windPow to be specified")
			( "windXi", po::value<double>(), "Ionization parameter, the ratio of the radiation and gas pressures" )
			( "windTic", po::value<double>(), "Inverse Compton temperature, K. Characterizes the hardness of the irradiating spectrum")
			( "windPow", po::value<double>(), "Multiplicative coefficient to control wind power")
			;
	return od;
}


SelfIrradiationOptions::SelfIrradiationOptions(const po::variables_map &vm, const DiskStructureArguments &dsa_args):
		SelfIrradiationArguments(
				vm["Cirr"].as<double>(),
				vm["irrindex"].as<double>(),
				vm["Cirrcold"].as<double>(),
				vm["irrindexcold"].as<double>(),
				vm["h2rcold"].as<double>(),
				vm["angulardistdisk"].as<std::string>()) {
	if (Cirr <= 0. && dsa_args.boundcond == "Tirr") {
		throw po::error("Set positive --Cirr when --boundcond=Tirr");
	}
}

po::options_description SelfIrradiationOptions::description() {
	po::options_description od("Parameters of self-irradiation.\nQirr = Cirr * (H/r / 0.05)^irrindex * L * psi / (4 pi R^2), where psi is angular distrbution of X-ray radiation");
	od.add_options()
			( "Cirr", po::value<double>()->default_value(default_Cirr), "Irradiation factor for the hot disk" )
			( "irrindex", po::value<double>()->default_value(default_irrindex), "Irradiation index for the hot disk" )
			( "Cirrcold", po::value<double>()->default_value(default_Cirr_cold), "Irradiation factor for the cold disk" )
			( "irrindexcold", po::value<double>()->default_value(default_irrindex_cold), "Irradiation index for the cold disk" )
			( "h2rcold", po::value<double>()->default_value(default_height_to_radius_cold), "Seme-height to radius ratio for the cold disk, it affects disk shadow in star" )
			( "angulardistdisk", po::value<std::string>()->default_value(default_angular_dist_disk), "Angular distribution of the disk X-ray radiation. Values: isotropic, plane" )
			;
	return od;
}


FluxOptions::FluxOptions(const po::variables_map &vm):
		FluxArguments(
				vm["colourfactor"].as<double>(),
				kevToHertz(vm["emin"].as<double>()),
				kevToHertz(vm["emax"].as<double>()),
				vm["staralbedo"].as<double>(),
				vm["inclination"].as<double>(),
				dayToS(vm["ephemerist0"].as<double>()),
				kpcToCm(vm["distance"].as<double>()),
				vm.count("colddiskflux") > 0,
				vm.count("starflux") > 0,
				lambdasInitializer(vm),
				passbandsInitializer(vm)) {}

vecd FluxOptions::lambdasInitializer(const po::variables_map &vm) {
	if (vm.count("lambda") == 0) {
		return vecd();
	}
	vecd lambdas(vm["lambda"].as<vecd>());
	transform(lambdas.begin(), lambdas.end(), lambdas.begin(), angstromToCm);
	return lambdas;
}

std::vector<Passband> FluxOptions::passbandsInitializer(const po::variables_map& vm) {
	if (vm.count("passband") == 0) {
		return {};
	}
	auto filepaths = vm["passband"].as<std::vector<std::string>>();
	std::vector<Passband> passbands;
	for (const auto &filepath : filepaths) {
		try {
			passbands.emplace_back(filepath);
		} catch (const std::ios_base::failure& e) {
			throw po::invalid_option_value("Passband file doesn't exist");
		}
	}
	return passbands;
}

po::options_description FluxOptions::description() {
	po::options_description od("Parameters of flux calculation");
	od.add_options()
			( "colourfactor", po::value<double>()->default_value(default_colourfactor), "Colour factor to calculate X-ray flux"  )
			( "emin", po::value<double>()->default_value(hertzToKev(default_emin)), "Minimum energy of X-ray band, keV" )
			( "emax", po::value<double>()->default_value(hertzToKev(default_emax)), "Maximum energy of X-ray band, keV" )
			( "staralbedo", po::value<double>()->default_value(default_star_albedo), "Part of X-ray radiation reflected by optical star, (1 - albedo) heats star's photosphere. Used only when --starflux is specified" )
			( "inclination,i", po::value<double>()->default_value(default_inclination), "Inclination of the system, degrees" )
			( "ephemerist0", po::value<double>()->default_value(default_ephemeris_t0), "Ephemeris for the time of the minimum of the orbital light curve T0, phase zero corresponds to inferior conjunction of the optical star, days" )
			( "distance", po::value<double>()->required(), "Distance to the system, kpc" )
			( "colddiskflux", "Add Fnu for cold disk into output file. Default output is for hot disk only" )
			( "starflux", "Add Fnu for irradiated optical star into output file. See --Topt, --starlod and --h2rcold options. Default is output for the hot disk only" )
			( "lambda", po::value<vecd>()->multitoken()->composing(), "Wavelength to calculate Fnu, Angstrom. You can use this option multiple times. For each lambda one additional column with values of spectral flux density Fnu [erg/s/cm^2/Hz] is produced" )
			( "passband", po::value<std::vector<std::string>>()->multitoken()->composing(), "Path of a file containing tabulated passband, the first column for wavelength in Angstrom, the second column for transmission factor, columns should be separated by spaces" )
			;
	return od;
}


CalculationOptions::CalculationOptions(const po::variables_map &vm):
		CalculationArguments(
				dayToS(vm["inittime"].as<double>()),
				dayToS(vm["time"].as<double>()),
				tauInitializer(vm),
				vm["Nx"].as<unsigned int>(),
				vm["gridscale"].as<std::string>(),
				vm["starlod"].as<unsigned int>()) {
	if (gridscale != "log" && gridscale != "linear") {
		throw po::invalid_option_value("Invalid --gridscale value");
	}
}

std::optional<double> CalculationOptions::tauInitializer(const po::variables_map& vm) {
	if (vm.count("tau")) {
		return dayToS(vm["tau"].as<double>());
	}
	return {};
}

po::options_description CalculationOptions::description() {
	po::options_description od("Parameters of disk evolution calculation");
	od.add_options()
			("inittime", po::value<double>()->default_value(default_init_time), "Initial time moment, days" )
			( "time,T", po::value<double>()->required(), "Time interval to calculate evolution, days" )
			( "tau",	po::value<double>(), "Time step, days" )
			( "Nx",	po::value<unsigned int>()->default_value(default_Nx), "Size of calculation grid" )
			( "gridscale", po::value<std::string>()->default_value(default_gridscale), "Type of grid for angular momentum h: log or linear" )
			( "starlod", po::value<unsigned int>()->default_value(default_starlod), "Level of detail of the optical star 3-D model. The optical star is represented by a triangular tile, the number of tiles is 20 * 4^starlod" )
			;
	return od;
}



FreddiOptions::FreddiOptions(const po::variables_map& vm) {
	if ((vm.count("starflux") > 0)
		&& (vm.count("Mx") == 0 || vm.count("Topt") == 0 || vm.count("Mopt") == 0 || vm.count("period") == 0)) {
		throw po::invalid_option_value("--starflux requires --Mx, --Mopt and --period to be specified");
	}

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
