#include <algorithm>
#include <boost/any.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "nonlinear_diffusion.hpp"
#include "spectrum.hpp"
#include "opacity_related.hpp"
#include "orbit.hpp"
#include "gsl_const_cgsm.h"


namespace po = boost::program_options;
namespace odeint = boost::numeric::odeint;
using namespace std;
using namespace std::placeholders;



int main(int ac, char *av[]){
	const double DAY = 86400.;
	const double Angstrem = 1e-8;
	const double keV = 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
	const double Jy = 1e-23;
	const double solar_radius = 6.955e10;
	const double kpc = 1000. * GSL_CONST_CGSM_PARSEC;

	// Allen's Astrophysical Quantities (4th ed.)
	const double lambdaU = 3600. * Angstrem;
	const double irr0U = 4.22e-9 / Angstrem;
	const double lambdaB = 4400. * Angstrem;
	const double irr0B = 6.4e-9 / Angstrem;
	const double lambdaV = 5500. * Angstrem;
	const double irr0V = 3.750e-9 / Angstrem;
	const double lambdaR = 7100 * Angstrem;
	const double irr0R = 1.75e-9 / Angstrem;
	const double lambdaI = 9700 * Angstrem;
	const double irr0I = 0.84e-9 / Angstrem;
	// Campins et al., 1985, AJ, 90, 896
	const double lambdaJ = 12600 * Angstrem;
	const double irr0J = 1600 * Jy *  GSL_CONST_CGSM_SPEED_OF_LIGHT / (lambdaJ*lambdaJ); 

	po::variables_map po_vm;

	double alpha = 0.25;
	double fc = 1.7;
	double kerr = 0.;
	double Mx = 5. * GSL_CONST_CGSM_SOLAR_MASS;
	double Mopt = 0.5 * GSL_CONST_CGSM_SOLAR_MASS;
	double P = 0.25 * DAY;
	double inclination = 0.;  // degrees
	double Distance = 10. * kpc;
	double r_in = 0.;
	double r_out;
	double T_min_hot_disk = 0.;
//	double k_irr = 0.05; //0.05; // (dlog H / dlog r - 1)
	double C_irr_input = 0.; // 1e-4;
	double mu = 0.62;
	double nu_min = 1. * keV;
	double nu_max = 12. * keV;
	int Nx = 1000;
	string grid_scale = "log";
	double Time = 50. * DAY;
	double tau = 0.25 * DAY;
	double eps = 1e-6;
	string bound_cond_type = "Teff";
	double F0 = 2e38;
	double sigma_for_F_gauss = 5.;
	double r_gauss_cut_to_r_out = 0.01;
	double gauss_sigma = 0.25;
	double gauss_mu = 1.;
	double power_order = 6.;
	double kMdot_out = 2.;
	string filename_prefix = "freddi";
	string output_dir = ".";
	bool output_fulldata = false;
	string initial_cond_shape = "powerF";
	string opacity_type = "Kramers";
	string irr_factor_type = "const";
	vector<double> additional_lambdas;

	const string config_filename = "freddi.ini";
	const char* home = getenv("HOME");
	const array<const string, 4> path_config_file = {".", home, INSTALLPATHPREFIX"/etc", "/etc"};

	double Mdot_in = 0.;
	double Mdot_in_prev;
	double Mdot_out = 0.;

	double Mdisk = 0.;

	{
		po::options_description desc("Freddi - numerical calculation of accretion disk evolution");

		po::options_description general("General options");
		general.add_options()
			( "help,h", "Produce help message" )
			( "prefix", po::value<string>(&filename_prefix)->default_value(filename_prefix), "Set prefix for output filenames. Output file with distribution of parameters over time is PREFIX.dat" )
			( "dir,d", po::value<string>(&output_dir)->default_value(output_dir), "Choose the directory to write output files. It should exist" )
			( "fulldata", "Output files PREFIX_%d.dat with radial structure for every time step. Default is to output only PREFIX.dat with global disk parameters for every time step" )
		;
		desc.add(general);
		
		po::options_description binary("Basic binary and disk parameters");
		binary.add_options()
			( "Mx,M", po::value<double>()->default_value(Mx/GSL_CONST_CGSM_SOLAR_MASS), "Mass of the central object, in the units of solar masses" )
			( "kerr", po::value<double>(&kerr)->default_value(kerr), "Dimensionless Kerr parameter of the black hole" )
			( "alpha,a", po::value<double>(&alpha)->default_value(alpha), "Alpha parameter of Shakura-Sunyaev model" )
			( "rin", po::value<double>(), "Inner radius of the disk, in the units of the Schwarzschild radius of the central object 2GM/c^2. If it isn't set then the radius of ISCO orbit is used defined by --Mx and --kerr values" )
			( "Mopt",	po::value<double>()->default_value(Mopt/GSL_CONST_CGSM_SOLAR_MASS), "Mass of the optical star, in units of solar masses" )
			( "period,P", po::value<double>()->default_value(P/DAY), "Orbital period of the binary system, in units of days" )
			( "rout,R", po::value<double>(), "Outer radius of the disk, in units of solar radius. If it isn't set then the tidal radius is used defined by --Mx, --Mopt and --period values" )
		;
		desc.add(binary);

		po::options_description internal("Parameters of the disk model");
		internal.add_options()
			( "opacity,O", po::value<string>(&opacity_type)->default_value(opacity_type), "Opacity law: Kramers (varkappa ~ rho / T^7/2) or OPAL (varkappa ~ rho / T^5/2)" )
			( "boundcond", po::value<string>(&bound_cond_type)->default_value(bound_cond_type), "Outer boundary movement condition\n\n"
				"Values:\n"
				"  Teff: outer radius of the disk moves inwards to keep photosphere temperature of the disk larger than some value. This value is specified by --Thot option\n"
				"  Tirr: outer radius of the disk moves inwards to keep irradiation flux of the disk larger than some value. The value of this minimal irradiation flux is [Stefan-Boltzmann constant] * Tirr^4, where Tirr is specified by --Thot option" ) // fourSigmaCrit, MdotOut
			( "Thot", po::value<double>(&T_min_hot_disk)->default_value(T_min_hot_disk), "Minimum photosphere or irradiation temperature at the outer edge of the hot disk, Kelvin. For details see --boundcond description" )
			( "F0", po::value<double>(&F0)->default_value(F0), "Initial maximum viscous torque in the disk, dyn*cm. Can be overwritten via --Mdisk0 and --Mdot0" )
			( "Mdisk0", po::value<double>(&Mdisk), "Initial disk mass, g. If both --F0 and --Mdisk0 are specified then --Mdisk0 is used. If both --Mdot0 and --Mdisk0 are specified then --Mdot0 is used" )
			( "Mdot0", po::value<double>(&Mdot_in), "Initial mass accretion rate through the inner radius, g/s. If --F0, --Mdisk0 and --Mdot0 are specified then --Mdot0 is used. Works only when --initialcond is set to sinusF or quasistat" )
			( "initialcond", po::value<string>(&initial_cond_shape)->default_value(initial_cond_shape), "Type of the initial condition for viscous torque F or surface density Sigma\n\n"
				"Values:\n"
				"  powerF: F ~ xi^powerorder, powerorder is specified by --powerorder option\n" // power option does the same
				"  powerSigma: Sigma ~ xi^powerorder, powerorder is specified by --powerorder option\n"
				"  sinusF: F ~ sin( xi * pi/2 )\n" // sinus option does the same
				"  gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2), mu and sigma are specified by --gaussmu and --gausssigma options\n"
				"  quasistat: F ~ f(h/h_out) * xi * h_out/h, where f is quasi-stationary solution found in Lipunova & Shakura 2000. f(xi=0) = 0, df/dxi(xi=1) = 0\n\n"
				"Here xi is (h - h_in) / (h_out - h_in)\n") // sinusparabola, sinusgauss
			( "powerorder", po::value<double>(&power_order)->default_value(power_order), "Parameter for the powerlaw initial condition distribution. This option works only with --initialcond=powerF or powerSigma" )
			( "gaussmu", po::value<double>(&gauss_mu)->default_value(gauss_mu), "Position of the maximum for Gauss distribution, positive number not greater than unity. This option works only with --initialcond=gaussF" )
			( "gausssigma", po::value<double>(&gauss_sigma)->default_value(gauss_sigma), "Width of for Gauss distribution. This option works only with --initialcond=gaussF" )
		;
		desc.add(internal);

		po::options_description irradiation("Parameters of self-irradiation");
		irradiation.add_options()
			( "Cirr", po::value<double>(&C_irr_input)->default_value(C_irr_input), "Irradiation factor" )
			( "irrfactortype", po::value<string>(&irr_factor_type)->default_value(irr_factor_type), "Type of irradiation factor Cirr\n\n"
				"Values:\n"
				"  const: doesn't depend on disk shape:\n[rad. flux] = Cirr  L / (4 pi r^2)\n"
				"  square: Cirr depends on the disk relative half-thickness:\n[rad. flux] = Cirr (z/r)^2 L / (4 pi r^2)\n\n"
				"Here L is bolometric Luminosity:\nL = eta Mdot c^2" )
		;
		desc.add(irradiation);

		po::options_description emission("Parameters of optical magnitudes calculation");
		emission.add_options()
			( "colourfactor", po::value<double>(&fc)->default_value(fc), "Colour factor to calculate X-ray flux"  )
			( "emin", po::value<double>()->default_value(nu_min/keV), "Minimum energy of X-ray band, keV" )
			( "emax", po::value<double>()->default_value(nu_max/keV), "Maximum energy of X-ray band, keV" )
			( "inclination,i", po::value<double>(&inclination)->default_value(inclination), "Inclination of the system, degrees" )
			( "distance", po::value<double>()->default_value(Distance/kpc), "Distance to the system, kpc" )
			( "lambda", po::value< vector<double> >(&additional_lambdas)->multitoken(), "Wavelength to calculate Fnu, Angstrom. You can use this option multiple times. For each lambda one additional column with values of spectral flux density Fnu [erg/s/cm^2/Hz] is produced" )
		;
		desc.add(emission);

		po::options_description numeric("Parameters of disk evolution calculation");
		numeric.add_options()
			( "time,T", po::value<double>()->default_value(Time/DAY), "Time interval to calculate evolution, days" )
			( "tau",	po::value<double>()->default_value(tau/DAY), "Time step, days" )
			( "Nx",	po::value<int>(&Nx)->default_value(Nx), "Size of calculation grid" )
			( "gridscale", po::value<string>(&grid_scale)->default_value(grid_scale), "Type of grid for angular momentum h: log or linear" )
		;
		desc.add(numeric);

		try {
			po::store( po::parse_command_line(ac, av, desc), po_vm );
			for (const auto &path : path_config_file){
				ifstream config(path + "/" + config_filename);
				po::store( po::parse_config_file(config, desc), po_vm );
			}
			po::notify(po_vm);
		} catch (exception &e){
			cerr << "Error: " << e.what() << endl;
			return 1;
		}

		if ( po_vm.count("help") ){
			cout << desc << endl;
			return 0;
		}
		output_fulldata = po_vm.count("fulldata");
		Mopt = po_vm["Mopt"].as<double>() * GSL_CONST_CGSM_SOLAR_MASS;
		Mx = po_vm["Mx"].as<double>() * GSL_CONST_CGSM_SOLAR_MASS;
		P = po_vm["period"].as<double>() * DAY;
		Distance = po_vm["distance"].as<double>() * kpc;
		nu_min = po_vm["emin"].as<double>() * keV;
		nu_max = po_vm["emax"].as<double>() * keV;
		tau = po_vm["tau"].as<double>() * DAY;
		Time = po_vm["time"].as<double>() * DAY;
		if ( po_vm.count("rout") ){
			r_out = po_vm["rout"].as<double>() * solar_radius;
		} else{
			r_out = r_out_func( Mx, Mopt, P );
			cout << "--rout isn't specified, using tidal radius " << r_out/solar_radius << " solar radii" << endl;
		}
		if ( po_vm.count("rin") ){
			r_in = po_vm["rin"].as<double>() * 3. * 2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
		} else{
			r_in = r_in_func( Mx, kerr );
		}
		if ( C_irr_input <= 0. and bound_cond_type == "Tirr" ){
			throw po::error("It is obvious to use nonpositive --Cirr with --boundcond=Tirr");
		}
		transform(additional_lambdas.begin(), additional_lambdas.end(), additional_lambdas.begin(), bind1st(multiplies<double>(), Angstrem));
	}

	const double GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx;
	const double eta = efficiency_of_accretion(kerr);
	const double h_in = sqrt( GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx  * r_in );
	const double h_out = sqrt( GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx  * r_out );
	const double cosi = cos( inclination / 180 * M_PI );
	const double cosiOverD2 = cosi / Distance / Distance;

	OpacityRelated *oprel;
	try{
		oprel = new OpacityRelated(opacity_type, Mx, alpha, mu);
	} catch (invalid_argument){
		throw po::invalid_option_value(opacity_type);
	}

	auto wunc = [oprel](
		const vector<double> &h, const vector<double> &F,
		int first, int last	
	) ->vector<double>{
		vector<double> W( first > 0 ? first : 0,  0. );
		for ( int i = first; i <= last; ++i ){
			W.push_back(
				pow(F.at(i), 1. - oprel->m) * pow(h.at(i), oprel->n) / (1. - oprel->m) / oprel->D
			);
		}
		return W;
	};

	// Equation from Lasota, Dubus, Kruk A&A 2008, Menou et al. 1999. Sigma_cr is from their fig 8 and connected to point where Mdot is minimal.
	auto Sigma_hot_disk = [alpha, Mx](double r) ->double{
		return 39.9 * pow(alpha/0.1, -0.80) * pow(r/1e10, 1.11) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -0.37);
	};

	vector<double> h(Nx), R(Nx);
	for ( int i = 0; i < Nx; ++i ){
		if ( grid_scale == "log" ){ 
			h.at(i) = h_in * pow( h_out/h_in, i/(Nx-1.) );
		} else if ( grid_scale == "linear" ){
			h.at(i) = h_in + (h_out - h_in) * i/(Nx-1.);
		} else{
			throw po::invalid_option_value(grid_scale);
		}
		R.at(i) = h.at(i) * h.at(i) / GM;
	}

	vector<double> F(Nx);
	if ( initial_cond_shape == "sinusgauss" ){
		const double F0_sinus = 1e-6 * F0;
		const double h_cut_for_F_gauss = h_out / sqrt(r_gauss_cut_to_r_out);
		const double F_gauss_cut = F0 * exp( - (h_cut_for_F_gauss-h_out)*(h_cut_for_F_gauss-h_out) / (2. * h_out*h_out/(sigma_for_F_gauss*sigma_for_F_gauss)) );
		for ( int i = 0; i < Nx; ++i ){
			double F_gauss = F0 * exp( - (h.at(i)-h_out)*(h.at(i)-h_out) / (2. * h_out*h_out/(sigma_for_F_gauss*sigma_for_F_gauss)) ) - F_gauss_cut;
			F_gauss = F_gauss >= 0 ? F_gauss : 0.;
			const double F_sinus =  F0_sinus * sin( (h.at(i) - h_in) / (h_out - h_in) * M_PI / 2. );
			F.at(i) = F_gauss + F_sinus;
		}
	} else if ( initial_cond_shape == "power" or initial_cond_shape == "powerF" ){
		if ( Mdot_in != 0. ){
			throw po::invalid_option_value("It is obvious to use --Mdot with --initialcond=powerF");
		}
		if ( Mdisk > 0. ){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = (1. - oprel->m) * power_order;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
				stepper,
				[a,b,x0]( const double &y, double &dydx, double x ){
					dydx = pow(x, a) * pow(x + x0, b);
				},
				integral, 0., 1., 0.01
			);
			F0 = pow( Mdisk * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral, 1. / (1. - oprel->m) );
		}
		for ( int i = 0; i < Nx; ++i ){
			F.at(i) = F0 * pow( (h.at(i) - h_in) / (h_out - h_in), power_order );
		}
	} else if ( initial_cond_shape == "powerSigma" ){
		if ( Mdot_in != 0. ){
			throw po::invalid_option_value("It is obvious to use --Mdot with --initialcond=powerSigma");
		}
		if ( Mdisk > 0. ){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = power_order;
			const double b = 3.;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
				stepper,
				[a,b,x0]( const double &y, double &dydx, double x ){
					dydx = pow(x, a) * pow(x + x0, b);
				},
				integral, 0., 1., 0.01
			);
			F0 = pow( Mdisk * (1. - oprel->m) * oprel->D * pow(h_out, 3. - oprel->n) / pow(h_out - h_in, 4.) / integral, 1. / (1. - oprel->m) );
		}
		for ( int i = 0; i < Nx; ++i ){
			const double Sigma_to_Sigmaout = pow( (h.at(i) - h_in) / (h_out - h_in), power_order );
			F.at(i) = F0 * pow( h.at(i) / h_out, (3. - oprel->n) / (1. - oprel->m) ) * pow( Sigma_to_Sigmaout, 1. / (1. - oprel->m) );
		}
	} else if ( initial_cond_shape == "sinus" or initial_cond_shape == "sinusF" ){
		if ( Mdot_in > 0. ){
			F0 = Mdot_in * (h_out - h_in) * 2./M_PI;
		} else if ( Mdisk > 0. ){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel->m;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
				stepper,
				[a,b,x0]( const double &y, double &dydx, double x ){
					dydx = pow(sin(x * M_PI_2), a) * pow(x + x0, b);
				},
				integral, 0., 1., 0.01
			);
			F0 = pow( Mdisk * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral, 1. / (1. - oprel->m) );
		}
		for ( int i = 0; i < Nx; ++i ){
			F.at(i) = F0 * sin( (h.at(i) - h_in) / (h_out - h_in) * M_PI_2 );
		}
	} else if ( initial_cond_shape == "sinusparabola" ){
		const double h_F0 = h_out * 0.9;
		const double delta_h = h_out - h_F0;

		F0 = 1.24e13 * pow(Sigma_hot_disk(R.at(Nx-1)), 10./7.) * pow(h.at(Nx-1), 22./7.) * pow(GM, -10./7.) * pow(alpha, 8./7.);

		Mdot_out = -kMdot_out * F0 / (h_F0 - h_in) * M_PI*M_PI;

		for ( int i = 0; i < Nx; ++i ){
			if ( h.at(i) < h_F0 ){
				F.at(i) = F0 * sin( (h.at(i) - h_in) / (h_F0 - h_in) * M_PI / 2. );
			} else{
				F.at(i) = F0 * ( 1. - kMdot_out / (h_F0-h_in) / delta_h * M_PI / 4. * (h.at(i) - h_F0)*(h.at(i) - h_F0) );
			}
		}
	} else if( initial_cond_shape == "quasistat" ){
		if ( Mdot_in > 0. ){
			F0 = Mdot_in * (h_out - h_in) / h_out * h_in / oprel->f_F(h_in/h_out);
		} else if ( Mdisk > 0. ){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double x0 = h_in / (h_out - h_in);
			const double x1 = h_in / h_out;
			double integral = 0.;
			integrate_adaptive(
				stepper,
				[x0,x1,oprel]( const double &y, double &dydx, double x ){
					dydx = pow(oprel->f_F(x * (1. - x1) + x1) * x / (x * (1. - x1) + x1), 1. - oprel->m) * pow(x + x0, oprel->n);
				},
				integral, 0., 1., 0.01
			);
			F0 = pow( Mdisk * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral, 1. / (1. - oprel->m) );
		}
		for ( int i = 0; i < Nx; ++i ){
			const double xi_LS2000 = h.at(i) / h_out;
			F.at(i) = F0 * oprel->f_F(xi_LS2000) * (1. - h_in / h.at(i)) / (1. - h_in / h_out);
		}
	} else if( initial_cond_shape == "gaussF" ){
		if ( gauss_mu <= 0. or gauss_mu > 1. ){
			throw po::invalid_option_value("--gaussmu value should be large than 0 and not large than 1");
		}
		if ( Mdot_in != 0. ){
			throw po::invalid_option_value("Usage of --Mdot with --initialcond=gaussF produces unstable results and it isn't motivated physically. Use --F0 or --Mdisk0 instead");
			//F0 = Mdot_in * (h_out - h_in) * gauss_sigma*gauss_sigma / gauss_mu * exp( gauss_mu*gauss_mu / (2. * gauss_sigma*gauss_sigma) );
		} else if ( Mdisk > 0. ){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel->m;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
				stepper,
				[a,b,x0,gauss_mu,gauss_sigma]( const double &y, double &dydx, double x ){
					dydx = exp( -(x - gauss_mu)*(x - gauss_mu) * a / (2. * gauss_sigma*gauss_sigma) ) * pow(x + x0, b);
				},
				integral, 0., 1., 0.01
			);
			F0 = pow( Mdisk * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral, 1. / (1. - oprel->m) );
		}
		for ( int i = 0; i < Nx; ++i ){
			const double xi = (h.at(i) - h_in) / (h_out - h_in);
			F.at(i) = F0 * exp( -(xi - gauss_mu)*(xi - gauss_mu) / (2. * gauss_sigma*gauss_sigma) );
		}
	} else{
		throw po::invalid_option_value(initial_cond_shape);
	}

	ofstream output_sum( output_dir + "/" + filename_prefix + ".dat" );
	output_sum << "#t    Mdot Mdisk Rhot Cirrout H2R   Teffout Tirrout Qiir2Qvisout Lx    mU  mB  mV  mR  mI  mJ ";
	for ( int i = 0; i < additional_lambdas.size(); ++i ){
		output_sum << " Fnu" << i;
		for ( double j = 0; j < 9 - log10(i+0.1); ++j ){
			output_sum << " ";
		}
	}
	output_sum << "\n";
	output_sum << "#days g/s  g     Rsun float   float K       K       float        erg/s mag mag mag mag mag mag";
	for ( int i = 0; i < additional_lambdas.size(); ++i ){
		output_sum << " erg/s/cm^2/Hz";
	}
	output_sum << "\n";
	for ( const auto &it : po_vm ){
		auto &value = it.second.value();
		if ( auto v = boost::any_cast<uint32_t>(&value) ){
			output_sum << "# "
			           << it.first.c_str()
			           << "="
			           << *v
					   << "\n";
		} else if ( auto v = boost::any_cast<string>(&value) ){
			output_sum << "# "
			           << it.first.c_str()
			           << "="
			           << *v
					   << "\n";
		} else if ( auto v = boost::any_cast<double>(&value) ){
			output_sum << "# "
			           << it.first.c_str()
			           << "="
			           << *v
					   << "\n";
		} else if ( auto v = boost::any_cast<int>(&value) ){
			output_sum << "# "
			           << it.first.c_str()
			           << "="
			           << *v
					   << "\n";
		} else if ( auto v = boost::any_cast< vector<double> >(&value) ){
			for ( int i = 0; i < v->size(); ++i ){
				output_sum << "# "
						   << it.first.c_str()
						   << "="
						   << v->at(i)
						   << "  # "
						   << i
						   << "\n";
				}
		} else {
			output_sum << "error\n";
			// throw po::invalid_option_value(it.first.c_str());
		}
	}
	output_sum << flush;
	if ( po_vm.count("rout") == 0 ){
		output_sum << "# --rout hadn't been specified, tidal radius " << r_out/solar_radius << " Rsun was used" << std::endl;
	}

	vector<double> W(Nx, 0.);
	W = wunc(h, F, 1, Nx-1);

	for( int i_t = 0; i_t <= Time/tau + 0.001; ++i_t ){
		const double t = i_t * tau;
		// cout << t/DAY << endl;

		vector<double> Tph(Nx, 0.), Tph_vis(Nx, 0.), Tph_X(Nx, 0.), Tirr(Nx,0.), Sigma(Nx, 0.), Height(Nx, 0.);

		Mdot_in_prev = Mdot_in;
		Mdot_in = ( F.at(1) - F.at(0) ) / ( h.at(1) - h.at(0) );

		double C_irr;
		for ( int i = 1; i < Nx; ++i ){
			Sigma.at(i) = W.at(i) * GM*GM / ( 4.*M_PI *  pow(h.at(i), 3.) );
			Height.at(i) = oprel->Height(R.at(i), F.at(i));
			Tph_vis.at(i) = GM * pow(h.at(i), -1.75) * pow( 3. / (8.*M_PI) * F.at(i) / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph_X.at(i) = fc * T_GR( R.at(i), kerr, Mx, Mdot_in, R.front() );

			double Qx;
			if ( irr_factor_type == "const" ){
				C_irr = C_irr_input;
				Qx = C_irr_input * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (4.*M_PI * R.at(i)*R.at(i));
			} else if ( irr_factor_type == "square" ){
				C_irr = C_irr_input * (Height.at(i) / R.at(i)) * (Height.at(i) / R.at(i));
				Qx = C_irr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (4.*M_PI * R.at(i)*R.at(i));
			} else{
				throw po::invalid_option_value(irr_factor_type);
			}
			Tirr.at(i) = pow( Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph.at(i) = pow( pow(Tph_vis.at(i), 4.) + Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
		}
		
		const double Lx = Luminosity( R, Tph_X, nu_min, nu_max, 100 ) / pow(fc, 4.);
	
		const double mU = -2.5 * log10( I_lambda(R, Tph, lambdaU) * cosiOverD2 / irr0U );
		const double mB = -2.5 * log10( I_lambda(R, Tph, lambdaB) * cosiOverD2 / irr0B );
		const double mV = -2.5 * log10( I_lambda(R, Tph, lambdaV) * cosiOverD2 / irr0V );
		const double mR = -2.5 * log10( I_lambda(R, Tph, lambdaR) * cosiOverD2 / irr0R );
		const double mI = -2.5 * log10( I_lambda(R, Tph, lambdaI) * cosiOverD2 / irr0I );
		const double mJ = -2.5 * log10( I_lambda(R, Tph, lambdaJ) * cosiOverD2 / irr0J );

		Mdisk = 0.;
		for ( int i = 0; i < Nx; ++i ){
			double stepR;
			if ( i == 0              ) stepR = R.at(i+1) - R.at(i  );
			if ( i == Nx-1           ) stepR = R.at(i  ) - R.at(i-1);
			if ( i > 1 and i < Nx-1  ) stepR = R.at(i+1) - R.at(i-1);
			Mdisk += 0.5 * Sigma.at(i) * 2.*M_PI * R.at(i) * stepR;
		}

		if (output_fulldata){
			ostringstream filename;
			filename << output_dir << "/" << filename_prefix << "_" << i_t << ".dat";
			ofstream output( filename.str() );
			output << "#h      R  F      Sigma  Teff Tvis Tirr Height" << "\n";
			output << "#cm^2/s cm dyn*cm g/cm^2 K    K    K    cm" << "\n";
			output << "# Time = " << t / DAY << " Mdot_in = " << Mdot_in << endl;
			for ( int i = 1; i < Nx; ++i ){
				output		<< h.at(i)
					<< "\t" << R.at(i)
					<< "\t" << F.at(i)
					<< "\t" << Sigma.at(i)
					<< "\t" << Tph.at(i)
					<< "\t" << Tph_vis.at(i)
					<< "\t" << Tirr.at(i)
					<< "\t" << Height.at(i)
					<< endl;
			}
		}

		output_sum		<< t / DAY
				<< "\t" << Mdot_in
				<< "\t" << Mdisk
				<< "\t" << R.at(Nx-1) / solar_radius
				<< "\t" << C_irr
				<< "\t" << Height.at(Nx-1) / R.at(Nx-1)
				<< "\t" << Tph.at(Nx-1)
				<< "\t" << Tirr.at(Nx-1)
				<< "\t" << pow( Tirr.at(Nx-1) / Tph_vis.at(Nx-1), 4. )
				<< "\t" << Lx
				<< "\t" << mU
				<< "\t" << mB
				<< "\t" << mV
				<< "\t" << mR
				<< "\t" << mI
				<< "\t" << mJ;
		for ( auto &lambda : additional_lambdas ){
			output_sum
			    << "\t" << I_lambda(R, Tph, lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * cosiOverD2;
		}
		output_sum      << endl;
		
				
		try{
			nonlenear_diffusion_nonuniform_1_2 (tau, eps, 0., Mdot_out, wunc, h, F);
			W = wunc(h, F, 1, Nx-1);
		} catch (runtime_error er){
			cout << er.what() << endl;
			break;
		}

		
		int ii = Nx;
		if (bound_cond_type == "MdotOut"){
			Mdot_out = - kMdot_out * Mdot_in;
			do{
				ii--;
			} while( Sigma.at(ii) < Sigma_hot_disk(R[ii]) );

		} else if (bound_cond_type == "fourSigmaCrit"){
			do{
				ii--;
				// Equation from Menou et al. 1999. Factor 4 is from their fig 8 and connected to point where Mdot = 0.
			} while( Sigma.at(ii) <  4 * Sigma_hot_disk(R[ii]) );
		} else if ( bound_cond_type == "Teff" ){
			do{
				ii--;
			} while( Tph.at(ii) < T_min_hot_disk );
		} else if (  bound_cond_type == "Tirr" ){
			if ( Mdot_in >= Mdot_in_prev  and ( initial_cond_shape == "power" or initial_cond_shape == "sinusgauss" ) ){
				do{
					ii--;
				} while( Tph.at(ii) < T_min_hot_disk );
			} else{
				do{
					ii--;
				} while( Tirr.at(ii) < T_min_hot_disk );
			}
		} else{
			throw po::invalid_option_value(bound_cond_type);
		}

		if ( ii < Nx-1 ){
			Nx = ii+1;
			// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
			h.resize(Nx);
		}
	}

	delete oprel;

	return 0;
}
