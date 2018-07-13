#include <algorithm>  // transform
#include <array>
#include <cstdlib>  // getenv
#include <exception>
#include <fstream>

#include <boost/numeric/odeint.hpp>

#include "arguments.hpp"


namespace odeint = boost::numeric::odeint;


constexpr const char GeneralArguments::default_prefix[];
constexpr const char GeneralArguments::default_dir[];

GeneralArguments::GeneralArguments(const po::variables_map& vm):
		prefix(vm["prefix"].as<std::string>()),
		dir(vm["dir"].as<std::string>()),
		fulldata(vm.count("fulldata") > 0) {
}

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


constexpr const double BasicDiskBinaryArguments::default_alpha;
constexpr const double BasicDiskBinaryArguments::default_Mx;
constexpr const double BasicDiskBinaryArguments::default_kerr;
constexpr const double BasicDiskBinaryArguments::default_Mopt;
constexpr const double BasicDiskBinaryArguments::default_period;

BasicDiskBinaryArguments::BasicDiskBinaryArguments(const po::variables_map &vm):
		alpha(vm["alpha"].as<double>()),
		Mx(sunToGram(vm["Mx"].as<double>())),
		kerr(vm["kerr"].as<double>()),
		Mopt(sunToGram(vm["Mopt"].as<double>())),
		period(dayToS(vm["period"].as<double>())),
		rin(rinInitializer(vm, Mx, kerr)),
		rout(routInitializer(vm, Mx, Mopt, period)) {
	if (rin >= rout) {
		throw po::invalid_option_value("rin should be smaller rout");
	}
}

double BasicDiskBinaryArguments::rinInitializer(const po::variables_map &vm, double Mx, double kerr) {
	if (vm.count("rin")) {
		return rgToCm(vm["rin"].as<double>(), Mx);
	}
	return rinFromMxKerr(Mx, kerr);
}

double BasicDiskBinaryArguments::routInitializer(const po::variables_map &vm, double Mx, double Mopt, double period) {
	if (vm.count("rout")) {
		return sunToCm(vm["rout"].as<double>());
	}
	return routFromMxMoptPeriod(Mx, Mopt, period);
}

po::options_description BasicDiskBinaryArguments::description() {
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
constexpr const double DiskStructureArguments::default_Thot;
constexpr const char DiskStructureArguments::default_boundcond[];
constexpr const double DiskStructureArguments::default_F0;
constexpr const double DiskStructureArguments::default_powerorder;
constexpr const double DiskStructureArguments::default_gaussmu;
constexpr const double DiskStructureArguments::default_gausssigma;
constexpr const double DiskStructureArguments::mu;

DiskStructureArguments::DiskStructureArguments(const po::variables_map &vm, const BasicDiskBinaryArguments& bdb_args):
		opacity(vm["opacity"].as<std::string>()),
		oprel(new OpacityRelated(opacity, bdb_args.Mx, bdb_args.alpha, mu)),
		boundcond(vm["boundcond"].as<std::string>()),
		Thot(vm["Thot"].as<double>()),
		initialcond(vm["initialcond"].as<std::string>()),
		Mdisk0(Mdisk0Initializer(vm)),
		Mdot0(Mdot0Initializer(vm)),
		powerorder(vm["powerorder"].as<double>()),
		gaussmu(vm["gaussmu"].as<double>()),
		gausssigma(vm["gausssigma"].as<double>()),
		F0(F0Initializer(vm, bdb_args)) {}

double DiskStructureArguments::Mdisk0Initializer(const po::variables_map& vm) {
	if (vm.count("Mdisk0") > 0) {
		return vm["Mdisk0"].as<double>();
	}
	return -1;
}

double DiskStructureArguments::Mdot0Initializer(const po::variables_map& vm) {
	if (vm.count("Mdot0") > 0) {
		return vm["Mdot0"].as<double>();
	}
	return -1;
}

double DiskStructureArguments::F0Initializer(const po::variables_map& vm, const BasicDiskBinaryArguments& bdb_args) {
	auto opacity = vm["opacity"].as<std::string>();
	OpacityRelated oprel(opacity, bdb_args.Mx, bdb_args.alpha, mu);
	auto initialcond = vm["initialcond"].as<std::string>();
	auto Mdisk0 = Mdisk0Initializer(vm);
	auto Mdot0 = Mdot0Initializer(vm);
	auto powerorder = vm["powerorder"].as<double>();
	auto gaussmu = vm["gaussmu"].as<double>();
	auto gausssigma = vm["gausssigma"].as<double>();
	auto F0 = vm["F0"].as<double>();

	const double h_in = bdb_args.h(bdb_args.rin);
	const double h_out = bdb_args.h(bdb_args.rout);

	const bool is_Mdot0_specified = (vm.count("Mdot0") > 0);
	const bool is_Mdisk0_specified = (vm.count("Mdisk0") > 0);

	if (initialcond == "powerF" || initialcond == "power") {
		if (is_Mdot0_specified) {
			throw po::invalid_option_value("It is obvious to use --Mdot with --initialcond=powerF");
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = (1. - oprel.m) * powerorder;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a,b,x0]( const double &y, double &dydx, double x ){
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0 = pow(Mdisk0 * (1. - oprel.m) * oprel.D / pow(h_out - h_in, oprel.n + 1.) / integral, 1. / (1. - oprel.m));
			return F0;
		}
		return F0;
	}
	if (initialcond == "powerSigma") {
		if (is_Mdot0_specified) {
			throw po::invalid_option_value("It is obvious to use --Mdot with --initialcond=powerSigma");
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = powerorder;
			const double b = 3.;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0 = pow(Mdisk0 * (1. - oprel.m) * oprel.D * pow(h_out, 3. - oprel.n) / pow(h_out - h_in, 4.) / integral, 1. / (1. - oprel.m));
			return F0;
		}
		return F0;
	}
	if (initialcond == "sinusF" || initialcond == "sinus") {
		if (is_Mdot0_specified) {
			F0 = Mdot0 * (h_out - h_in) * 2./M_PI;
			return F0;
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel.m;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(sin(x * M_PI_2), a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0 = pow(Mdisk0 * (1. - oprel.m) * oprel.D / pow(h_out - h_in, oprel.n + 1.) / integral, 1. / (1. - oprel.m));
			return F0;
		}
		return F0;
	}
	if (initialcond == "gaussF") {
		if (gaussmu <= 0. or gaussmu > 1.){
			throw po::invalid_option_value("--gaussmu value should be large than 0 and not large than 1");
		}
		if (is_Mdot0_specified) {
			throw po::invalid_option_value("Usage of --Mdot with --initialcond=gaussF produces unstable results and it isn't motivated physically. Use --F0 or --Mdisk0 instead");
			//F0 = Mdot_in * (h_out - h_in) * gauss_sigma*gauss_sigma / gauss_mu * exp( gauss_mu*gauss_mu / (2. * gauss_sigma*gauss_sigma) );
		}
		if (is_Mdisk0_specified){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel.m;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a,b,x0,gaussmu,gausssigma]( const double &y, double &dydx, double x ){
						dydx = exp( -(x - gaussmu)*(x - gaussmu) * a / (2. * gausssigma*gausssigma) ) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0 = pow(Mdisk0 * (1. - oprel.m) * oprel.D / pow(h_out - h_in, oprel.n + 1.) / integral, 1. / (1. - oprel.m));
			return F0;
		}
		return F0;
	}
	if (initialcond == "quasistat") {
		if (is_Mdot0_specified) {
			F0 = Mdot0 * (h_out - h_in) / h_out * h_in / oprel.f_F(h_in/h_out);
			return F0;
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double x0 = h_in / (h_out - h_in);
			const double x1 = h_in / h_out;
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[x0,x1,&oprel]( const double &y, double &dydx, double x ){
						dydx = pow(oprel.f_F(x * (1. - x1) + x1) * x / (x * (1. - x1) + x1), 1. - oprel.m) * pow(x + x0, oprel.n);
					},
					integral, 0., 1., 0.01
			);
			F0 = pow(Mdisk0 * (1. - oprel.m) * oprel.D / pow(h_out - h_in, oprel.n + 1.) / integral, 1. / (1. - oprel.m));
			return F0;
		}
		return F0;
	}
	throw po::invalid_option_value(initialcond);
}

po::options_description DiskStructureArguments::description() {
	po::options_description od("Parameters of the disk mode");
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


constexpr const double SelfIrradiationArguments::default_Cirr;
constexpr const char SelfIrradiationArguments::default_irrfactortype[];

SelfIrradiationArguments::SelfIrradiationArguments(const po::variables_map &vm, const DiskStructureArguments &dsa_args):
		Cirr(vm["Cirr"].as<double>()),
		irrfactortype(vm["irrfactortype"].as<std::string>()) {
	if (Cirr <= 0. && dsa_args.boundcond == "Tirr") {
		throw po::error("It is obvious to use nonpositive --Cirr with --boundcond=Tirr");
	}
	if (irrfactortype != "const" && irrfactortype != "square") {
		throw po::invalid_option_value("--irrfactortype has invalid value");
	}
}

po::options_description SelfIrradiationArguments::description() {
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


constexpr const double FluxArguments::default_colourfactor;
constexpr const double FluxArguments::default_emin;
constexpr const double FluxArguments::default_emax;
constexpr const double FluxArguments::default_inclination;
constexpr const double FluxArguments::default_distance;

FluxArguments::FluxArguments(const po::variables_map &vm):
		colourfactor(vm["colourfactor"].as<double>()),
		emin(kevToHertz(vm["emin"].as<double>())),
		emax(kevToHertz(vm["emax"].as<double>())),
		inclination(vm["inclination"].as<double>()),
		distance(kpcToCm(vm["distance"].as<double>())),
		lambdas(lambdasInitializer(vm)) {}

std::vector<double> FluxArguments::lambdasInitializer(const po::variables_map &vm) const {
	if (vm.count("lambda") > 0) {
		std::vector<double> lambdas(vm["lambda"].as<std::vector<double> >());
		transform(lambdas.begin(), lambdas.end(), lambdas.begin(), angstromToCm);
		return lambdas;
	}
	return std::vector<double>();
}

po::options_description FluxArguments::description() {
	po::options_description od("Parameters of flux calculation");
	od.add_options()
			( "colourfactor", po::value<double>()->default_value(default_colourfactor), "Colour factor to calculate X-ray flux"  )
			( "emin", po::value<double>()->default_value(hertzToKev(default_emin)), "Minimum energy of X-ray band, keV" )
			( "emax", po::value<double>()->default_value(hertzToKev(default_emax)), "Maximum energy of X-ray band, keV" )
			( "inclination,i", po::value<double>()->default_value(default_inclination), "Inclination of the system, degrees" )
			( "distance", po::value<double>()->default_value(cmToKpc(default_distance)), "Distance to the system, kpc" )
			( "lambda", po::value< std::vector<double> >()->multitoken(), "Wavelength to calculate Fnu, Angstrom. You can use this option multiple times. For each lambda one additional column with values of spectral flux density Fnu [erg/s/cm^2/Hz] is produced" )
			;
	return od;
}


constexpr const double CalculationArguments::default_time;
constexpr const double CalculationArguments::default_tau;
constexpr const unsigned int CalculationArguments::default_Nx;
constexpr const char CalculationArguments::default_gridscale[];

CalculationArguments::CalculationArguments(const po::variables_map &vm):
		time(dayToS(vm["time"].as<double>())),
		tau(dayToS(vm["tau"].as<double>())),
		Nx(vm["Nx"].as<unsigned int>()),
		gridscale(vm["gridscale"].as<std::string>()),
		eps(1e-6) {
	if (gridscale != "log" && gridscale != "linear") {
		throw po::invalid_option_value("Invalid --gridscale value");
	}
}

po::options_description CalculationArguments::description() {
	po::options_description od("Parameters of disk evolution calculation");
	od.add_options()
			( "time,T", po::value<double>()->default_value(sToDay(default_time)), "Time interval to calculate evolution, days" )
			( "tau",	po::value<double>()->default_value(sToDay(default_tau)), "Time step, days" )
			( "Nx",	po::value<unsigned int>()->default_value(default_Nx), "Size of calculation grid" )
			( "gridscale", po::value<std::string>()->default_value(default_gridscale), "Type of grid for angular momentum h: log or linear" )
			;
	return od;
}


FreddiArguments::FreddiArguments(const po::variables_map& vm):
		general(new GeneralArguments(vm)),
		basic(new BasicDiskBinaryArguments(vm)),
		disk(new DiskStructureArguments(vm, *basic)),
		irr(new SelfIrradiationArguments(vm, *disk)),
		flux(new FluxArguments(vm)),
		calc(new CalculationArguments(vm)) {
}

po::options_description FreddiArguments::description() {
	po::options_description desc("Freddi - numerical calculation of accretion disk evolutio");
	desc.add(GeneralArguments::description());
	desc.add(BasicDiskBinaryArguments::description());
	desc.add(DiskStructureArguments::description());
	desc.add(SelfIrradiationArguments::description());
	desc.add(FluxArguments::description());
	desc.add(CalculationArguments::description());
	return desc;
}


po::variables_map parseArguments(int ac, char* av[]) {
	const std::string config_filename = "freddi.ini";
	const char* home = getenv("HOME");
	const std::array<const std::string, 4> path_config_file = {".", home, INSTALLPATHPREFIX"/etc", "/etc"};
	auto desc = FreddiArguments::description();
	po::variables_map vm;
	po::store( po::parse_command_line(ac, av, desc), vm );
	for (const auto &path : path_config_file){
		std::ifstream config(path + "/" + config_filename);
		po::store( po::parse_config_file(config, desc), vm );
	}
	po::notify(vm);
	return vm;
}
