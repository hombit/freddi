#include "arguments.hpp"


constexpr const char GeneralArguments::default_prefix[];
constexpr const char GeneralArguments::default_dir[];

po::options_description GeneralArguments::description() {
	po::options_description od("General options");
	od.add_options()
			( "help,h", "Produce help message" )
			( "prefix", po::value<std::string>()->default_value(default_prefix), "Set prefix for output filenames. Output file with distribution of parameters over time is PREFIX.dat" )
			( "dir,d", po::value<std::string>()->default_value(default_dir), "Choose the directory to write output files. It should exist" )
			( "fulldata", "Output files PREFIX_%d.dat with radial structure for every time step. Default is to output only PREFIX.dat with global disk parameters for every time step" )
			;
	return od;
}


po::options_description BasicDiskBinaryArguments::description() {
	po::options_description od("General options");
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


double BlackHoleFunctions::rISCORg(const double kerr) {  // From «Black Hole Accretion Disks», A.44 (p. 530)
	const double Z1 = 1. + std::cbrt((1. - kerr * kerr)) * (std::cbrt((1. + kerr)) + std::cbrt((1. - kerr)));
	const double Z2 = std::sqrt(3. * kerr * kerr + Z1 * Z1);
	return 3. + Z2 - std::sqrt((3. - Z1) * (3. + Z1 + 2. * Z2));
}


double BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(const double MxToMopt) { // Eggleton, P. P. 1983, ApJ, 268, 368
	const double q = std::cbrt(MxToMopt);
	return 0.49 * q * q / (0.6 * q * q + std::log(1. + q));
}


constexpr const char DiskStructureArguments::default_opacity[];
constexpr const char DiskStructureArguments::default_initialcond[];
constexpr const char DiskStructureArguments::default_boundcond[];

po::options_description DiskStructureArguments::description() {
	po::options_description od("General options");
	od.add_options()
			( "opacity,O", po::value<std::string>()->default_value(default_opacity), "Opacity law: Kramers (varkappa ~ rho / T^7/2) or OPAL (varkappa ~ rho / T^5/2)" )
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


constexpr const char SelfIrradiationArguments::irrfactortype[];

po::options_description SelfIrradiationArguments::description() {
	po::options_description od("General options");
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


constexpr const char CalculationArguments::default_gridscale[];


FreddiArguments::FreddiArguments(int argc, const char *argv[]) {

}