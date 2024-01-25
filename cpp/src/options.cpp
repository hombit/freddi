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
	po::options_description od("General options:");
	od.add_options()
			( "help,h", "Produce help message\n" )
			( "config", po::value<std::string>(), "Set filepath for additional configuration file. There is no need to declare a configuration file with the default name freddi.ini\n" )
			( "prefix", po::value<std::string>()->default_value(default_prefix), "Set prefix for output filenames. Output file with distribution of parameters over time is PREFIX.dat\n" )
			( "stdout", "Output temporal distribution to stdout instead of PREFIX.dat file\n" )
			( "dir,d", po::value<std::string>()->default_value(default_dir), "Choose the directory to write output files. It should exist\n" )
			( "precision", po::value<unsigned int>()->default_value(default_output_precision), "Number of digits to print into output files\n" )
			( "tempsparsity", po::value<unsigned int>()->default_value(default_temp_sparsity_output), "Output every k-th time moment\n" )
			( "fulldata", "Output files PREFIX_%d.dat with radial structure for every time step. Default is to output only PREFIX.dat with global disk parameters for every time step\n" )
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
	po::options_description od("Basic binary and disk parameters\n");
	od.add_options()
			( "alpha,a", po::value<double>()->required(), "Alpha parameter of Shakura-Sunyaev model\n" )
			( "alphacold", po::value<double>(), "Alpha parameter of cold disk, currently it is used only for the critical maximum value of the surface density of the cold disk Sigma_minus (Lasota et al., 2008, A&A 486, 523) and the cooling front velocity (Ludwig et al., 1994, A&A 290, 473), see --Qirr2Qvishot. Default value is --alpha divided by ten\n" )  // default_alpha_to_alphacold
			( "Mx,M", po::value<double>()->required(), "Mass of the central object, in the units of solar masses\n" )
			( "kerr", po::value<double>()->default_value(default_kerr), "Dimensionless Kerr parameter of the black hole\n" )
			( "Mopt",	po::value<double>()->required(), "Mass of the optical star, in units of solar masses\n" )
			( "rochelobefill", po::value<double>()->default_value(default_roche_lobe_fill), "Dimensionless factor describing a size of the optical star. Polar radius of the star is rochelobefill * (polar radius of critical Roche lobe)\n" )
			( "Topt", po::value<double>()->default_value(default_Topt), "Effective temperature of the optical star, in units of Kelvins\n" )
			( "period,P", po::value<double>()->required(), "Orbital period of the binary system, in units of days\n" )
			( "rin", po::value<double>(), "Inner radius of the disk, in the units of the gravitational radius of the central object GM/c^2. There is no need to set it for a neutron star. If it isn't set for a black hole then the radius of ISCO orbit is used, defined by --Mx and --kerr values\n" )
			( "rout,R", po::value<double>(), "Outer radius of the disk, in units of solar radius. If it isn't set then the tidal radius is used, defined by --Mx, --Mopt and --period values as 90% of the Roche lobe radius (Papaloizou & Pringle, 1977, MNRAS, 181, 441; see also Artymowicz & Lubow, 1994, ApJ, 421, 651; http://xray.sai.msu.ru/~galja/images/tidal_radius.pdf)\n" )
			( "risco", po::value<double>(), "Innermost stable circular orbit, in units of gravitational radius of the central object GM/c^2. If it isn't set then the radius of ISCO orbit is used defined by --Mx and --kerr values\n" )
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
				vm["Rhot_Mdotzero_factor"].as<double>(),     
				vm["check_state_approach"].as<std::string>(),       
				vm["check_Sigma_approach"].as<std::string>(),
				vm["check_Temp_approach"].as<std::string>(),
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
	if (windtype == "ShieldsOscil1986"){
		if (vm.count("windC_w") == 0) {
			throw po::error("--windC_w is required if --windtype=ShieldsOscil1986");
		}
		if (vm.count("windR_w") == 0) {
			throw po::error("--windR_w is required if --windtype=ShieldsOscil1986");
		}
		return {
				{"C_w", vm["windC_w"].as<double>()},
				{"R_w", vm["windR_w"].as<double>()}
		};
	}
	if (windtype == "Janiuk2015"){
		if (vm.count("windA_0") == 0) {
			throw po::error("--windA_0 is required if --windtype=Janiuk2015");
		}
		if (vm.count("windB_1") == 0) {
			throw po::error("--windB_1 is required if --windtype=Janiuk2015");
		}
		return {
				{"A_0", vm["windA_0"].as<double>()},
				{"B_1", vm["windB_1"].as<double>()}
		};
	}
	if (windtype == "Shields1986"){
		if (vm.count("windXi_max") == 0) {
			throw po::error("--windXi_max is required if --windtype=Shields1986");
		}
		if (vm.count("windT_ic") == 0) {
			throw po::error("--windT_ic is required if --windtype=Shields1986");
		}
		if (vm.count("windPow") == 0) {
			throw po::error("--windPow is required if --windtype=Shields1986");
		}
		return {
				{"Xi_max", vm["windXi_max"].as<double>()},
				{"T_ic", vm["windT_ic"].as<double>()},
				{"Pow", vm["windPow"].as<double>()}
		};
	}
	if (windtype == "Woods1996AGN"){
		if (vm.count("windC_0") == 0) {
			throw po::error("--windC_0 is required if --windtype=Woods1996AGN");
		}
		if (vm.count("windT_ic") == 0) {
			throw po::error("--windT_ic is required if --windtype=Woods1996AGN");
		}
		return {
				{"C_0", vm["windC_0"].as<double>()},
				{"T_ic", vm["windT_ic"].as<double>()}
		};
	}
	if (windtype == "Woods1996"){
		if (vm.count("windXi_max") == 0) {
			throw po::error("--windXi_max is required if --windtype=Woods1996");
		}
		if (vm.count("windT_ic") == 0) {
			throw po::error("--windT_ic is required if --windtype=Woods1996");
		}
		if (vm.count("windPow") == 0) {
			throw po::error("--windPow is required if --windtype=Woods1996");
		}
		return {
				{"Xi_max", vm["windXi_max"].as<double>()},
				{"T_ic", vm["windT_ic"].as<double>()},
				{"Pow", vm["windPow"].as<double>()}
		};
	}
	if (windtype == "toy"){
		if (vm.count("windC_w") == 0) {
			throw po::error("--windC_w is required if --windtype=toy");
		}
		return {

				{"C_w", vm["windC_w"].as<double>()}
		};
	}
	

	throw po::invalid_option_value("Unknown --windtype=" + windtype);
}


po::options_description DiskStructureOptions::description() {
	po::options_description od(DiskStructureOptions::caption);
	od.add_options()
			( "opacity,O", po::value<std::string>()->default_value(default_opacity), "Opacity law: Kramers (varkappa ~ rho / T^7/2) or OPAL (varkappa ~ rho / T^5/2)\n" )
			( "Mdotout", po::value<double>()->default_value(default_Mdotout), "Accretion rate onto the disk through its outer radius\n" )
			( "boundcond", po::value<std::string>()->default_value(default_boundcond),
					"Outer-boundary movement condition\n\n"
					"Values:\n"
					"  Teff: outer radius of the disk moves inwards to keep photosphere temperature of the disk larger than some value. This value is specified by --Thot option\n"
					"  Tirr: outer radius of the disk moves inwards to keep irradiation flux of the disk larger than some value. The value of this minimal irradiation flux is [Stefan-Boltzmann constant] * Tirr^4, where Tirr is specified by --Thot option\n" ) // fourSigmaCrit, MdotOut
			( "Thot", po::value<double>()->default_value(default_Thot), "Minimum photosphere or irradiation temperature at the outer edge of the hot disk, Kelvin. For details see --boundcond description\n" )
			( "Qirr2Qvishot", po::value<double>()->default_value(m::pow<4>(default_Tirr2Tvishot)), "Minimum Qirr / Qvis ratio at the outer edge of the hot disk to switch the control over the evolution of the hot disk radius: from temperature-based regime to Sigma-based cooling-front regime (see Lipunova et al. (2021, Section 2.4) and Eq. A.1 in Lasota et al. 2008; --alpha value is used for Sigma_plus and --alphacold value is used for Sigma_minus)\n" )
			("Rhot_Mdotzero_factor", po::value<double>()->default_value(default_Rhot_Mdotzero_factor), "We check conditions for cooling front at current radius mpltiplied by Rhot_Mdotzero_factor\n" )
			("check_state_approach", po::value<std::string>()->default_value(default_check_state_approach), "Type of checking whether the ring is hot or cold\n\n"
					"Values:\n"
					" before2024: original version, as published in Lipunova&Malanchev (2017); Lipunova et al (2022); Avakyan et al (2024)\n"
					" logic: included option for checking conditions at radius different from the radius where accretion rate is zero. See Rhot_Mdotzero_factor check_Sigma_approach, and check_Temp_approach\n")
			("check_Sigma_approach", po::value<std::string>()->default_value(default_check_Sigma_approach), "Type of checking Sigma for hot or cold state\n\n"
					"Values:\n"
					" simple: assume that Sigma is proportional to R^(-3/4) between radius where Mdot = 0 and the cooling fron radius\n"   
					" Menou99a: assume that Sigma is 4.3 times less at the cooling front comparing to radius where Mdot = 0; See fig.8 of Menou et al. (1999 MNRAS 305, 79)\n" )
			("check_Temp_approach", po::value<std::string>()->default_value(default_check_Temp_approach), "Type of checking irradiation temperature for hot or cold state\n\n"
					"Values:\n"
					" const: assume that critical Tirr is constant, see --Thot, and --boundcond\n"   
					" Tavleed: assume that critical Tirr depends on Qvis/Qirr, see Tavleev et al (2023)\n" )
			( "initialcond", po::value<std::string>()->default_value(default_initialcond),
					"Type of the initial condition for viscous torque F or surface density Sigma\n\n"
					"Values:\n"
					"  [NB! Here below dimensionless xi = (h - h_in) / (h_out - h_in)]\n\n"
					"  powerF: F ~ xi^powerorder, powerorder is specified by --powerorder option\n" // power does the same
					"  linearF: F ~ xi, specific case of powerF but can be normalised by --Mdot0, see its description for details\n" // linear does the same
					"  powerSigma: Sigma ~ xi^powerorder, powerorder is specified by --powerorder option\n"
					"  sineF: F ~ sin( xi * pi/2 )\n" // sinus option does the same
					"  gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2), mu and sigma are specified by --gaussmu and --gausssigma options\n"
					"  quasistat: F ~ f(h/h_out) * xi * h_out/h, where f is quasi-stationary solution found in Lipunova & Shakura 2000. f(xi=0) = 0, df/dxi(xi=1) = 0\n")
			( "F0", po::value<double>(), "Initial maximum viscous torque in the disk, dyn*cm. Can be overwritten via --Mdisk0 and --Mdot0\n" )
			( "Mdisk0", po::value<double>(), "Initial disk mass, g. If both --F0 and --Mdisk0 are specified then --Mdisk0 is used. If both --Mdot0 and --Mdisk0 are specified then --Mdot0 is used\n" )
			( "Mdot0", po::value<double>(), "Initial mass accretion rate through the inner radius, g/s. If --F0, --Mdisk0 and --Mdot0 are specified then --Mdot0 is used. Works only when --initialcond is set to linearF, sinusF or quasistat\n" )
			( "powerorder", po::value<double>(), "Parameter for the powerlaw initial condition distribution. This option works only with --initialcond=powerF or powerSigma\n" )
			( "gaussmu", po::value<double>(), "Position of the maximum for Gauss distribution, positive number not greater than unity. This option works only with --initialcond=gaussF\n" )
			( "gausssigma", po::value<double>(), "Width of for Gauss distribution. This option works only with --initialcond=gaussF\n" )
			( "windtype", po::value<std::string>()->default_value(default_wind),
			        "Type of the wind\n\n"
					"  no: no wind\n"
					"  SS73C: super-Eddington spherical wind from Shakura-Sunyaev 1973\n"
					"  ShieldsOscil1986: toy wind model from Shields et al. 1986 which was used to obtain oscillations in the disk luminosity. Requires --windC_w and --windR_w to be specified\n"
					"  Janiuk2015: super-Eddington wind from Janiuk et al 2015. Requires --windA_0 and --windB_1 to be specified\n"
					"  Shields1986: thermal wind from Begelman et al. 1983 and Shields et al. 1986. Requires --windXi_max, --windT_ic and --windPow to be specified\n"
					"  Woods1996AGN: thermal AGN wind from Woods et al. 1996. Requires --windC_0 and --windT_ic to be specified\n"
					"  Woods1996: thermal wind from Woods et al. 1996. Requires --windXi_max, --windT_ic and --windPow to be specified\n"
					"  toy: a toy wind model used in arXiv:2105.11974, the mass loss rate is proportional to the central accretion rate. Requires --windC_w to be specified\n")
			( "windC_w", po::value<double>(), "The ratio of the mass loss rate due to wind to the central accretion rate, |Mwind|/Macc\n")
			( "windR_w", po::value<double>(), "The ratio of the wind launch radius to the outer disk radius, Rwind/Rout\n")		
			( "windA_0", po::value<double>(), "Dimensionless parameter characterizing the strength of the super-Eddington wind in the framework of the model Janiuk et al. 2015. Effective value range from 10 to 25\n")
			( "windB_1", po::value<double>(), "The quantity is of the order of unity. Characterizes the relationship between the change in energy per particle and virial energy.\nE = B_1 * k * T\n")
			( "windXi_max", po::value<double>(), "Ionization parameter, the ratio of the radiation and gas pressures\n" )
			( "windT_ic", po::value<double>(), "Inverse Compton temperature, K. Characterizes the hardness of the irradiating spectrum\n")
			( "windPow", po::value<double>(), "Multiplicative coefficient to control wind power\n")
			( "windC_0", po::value<double>(), "Characteristic column density of the wind mass loss rate from Woods et al. 1996 model, g/(s*cm^2). For AGN approx value is 3e-13 g/(s*cm^2)\n")
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
	po::options_description od("Parameters of self-irradiation:\nQirr = Cirr * (H/r / 0.05)^irrindex * L * psi / (4 pi R^2), where psi is the angular distribution of X-ray radiation\n");
	od.add_options()
			( "Cirr", po::value<double>()->default_value(default_Cirr), "Irradiation factor for the hot disk\n" )
			( "irrindex", po::value<double>()->default_value(default_irrindex), "Irradiation index for the hot disk\n" )
			( "Cirrcold", po::value<double>()->default_value(default_Cirr_cold), "Irradiation factor for the cold disk\n" )
			( "irrindexcold", po::value<double>()->default_value(default_irrindex_cold), "Irradiation index for the cold disk\n" )
			( "h2rcold", po::value<double>()->default_value(default_height_to_radius_cold), "Semi-height to radius ratio for the cold disk\n" )
			( "angulardistdisk", po::value<std::string>()->default_value(default_angular_dist_disk), "Angular distribution of the disk X-ray radiation. Values: isotropic, plane\n" )
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
	po::options_description od("Parameters of flux calculation:\n");
	od.add_options()
			( "colourfactor", po::value<double>()->default_value(default_colourfactor), "Colour factor to calculate X-ray flux\n"  )
			( "emin", po::value<double>()->default_value(hertzToKev(default_emin)), "Minimum energy of X-ray band, keV\n" )
			( "emax", po::value<double>()->default_value(hertzToKev(default_emax)), "Maximum energy of X-ray band, keV\n" )
			( "staralbedo", po::value<double>()->default_value(default_star_albedo), "Part of X-ray radiation reflected by optical star, (1 - albedo) heats star's photosphere. Used only when --starflux is specified\n" )
			( "inclination,i", po::value<double>()->default_value(default_inclination), "Inclination of the system, degrees\n" )
			( "ephemerist0", po::value<double>()->default_value(default_ephemeris_t0), "Ephemeris for the time of the minimum of the orbital light curve T0, phase zero corresponds to inferior conjunction of the optical star, days\n" )
			( "distance", po::value<double>()->required(), "Distance to the system, kpc\n" )
			( "colddiskflux", "Add Fnu for cold disk into output file. Default output is for hot disk only\n" )
			( "starflux", "Add Fnu for irradiated optical star into output file. See --Topt, --starlod and --h2rcold options. Default is output for the hot disk only\n" )
			( "lambda", po::value<vecd>()->multitoken()->composing(), "Wavelength to calculate Fnu, Angstrom. You can use this option multiple times. For each lambda one additional column with values of spectral flux density Fnu [erg/s/cm^2/Hz] is produced\n" )
			( "passband", po::value<std::vector<std::string>>()->multitoken()->composing(), "Path of a file containing tabulated passband, the first column for wavelength in Angstrom, the second column for transmission factor, columns should be separated by spaces\n" )
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
	po::options_description od("Parameters of disk evolution calculation:\n");
	od.add_options()
			("inittime", po::value<double>()->default_value(default_init_time), "Initial time moment, days\n" )
			( "time,T", po::value<double>()->required(), "Time interval to calculate evolution, days\n" )
			( "tau",	po::value<double>(), "Time step, days\n" )
			( "Nx",	po::value<unsigned int>()->default_value(default_Nx), "Size of calculation grid\n" )
			( "gridscale", po::value<std::string>()->default_value(default_gridscale), "Type of grid for angular momentum h: log or linear\n" )
			( "starlod", po::value<unsigned int>()->default_value(default_starlod), "Level of detail of the optical star 3-D model. The optical star is represented by a triangular tile, the number of tiles is 20 * 4^starlod\n" )
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
