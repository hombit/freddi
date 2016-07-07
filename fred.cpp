#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <gsl/gsl_const_cgsm.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "nonlinear_diffusion.hpp"
#include "spectrum.hpp"
#include "orbit.hpp"


namespace po = boost::program_options;
using namespace std;
using namespace std::placeholders;



int main(int ac, char *av[]){
	const double DAY = 86400.;
	const double Angstrem = 1e-8;
	const double keV = 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
	const double solar_radius = 6.955e10;
	const double kpc = 1000. * GSL_CONST_CGSM_PARSEC;

	const double lambdaB = 4400. * Angstrem;
	const double irr0B = 6.61e-9 / Angstrem;
	// Calculated for BUXTON & BAILYN 2004
	const double lambdaV = 5450. * Angstrem;
	const double irr0V = 2.90e-9 / Angstrem;

	double alpha = 0.55;
	double fc = 1.7;
	double kerr = 0.;
	double eta = efficiency_of_accretion(kerr);
	double Mx = 7.5 * GSL_CONST_CGSM_SOLAR_MASS;
	double Mopt = 0.8 * GSL_CONST_CGSM_SOLAR_MASS;
	double P = 0.433 * DAY;
	double cosi = cos(20.7 / 180 * M_PI);
	double Distance = 10. * kpc;
	double r_in = r_in_func( Mx, kerr );
	double r_out = r_out_func( Mx, Mopt, P );
	double T_min_hot_disk = 8000;
//	double k_irr = 0.05; //0.05; // (dlog H / dlog r - 1)
	double C_irr = 0.; // 1e-4;
	double mu = 0.5;
	double nu_min = 1.2 * keV;
	double nu_max = 37.2 * keV;
	int Nx = 1000;
	string grid_scale = "log";
	double Time = 30. * DAY;
	double tau = 0.25 * DAY;
	double eps = 1e-6;
	string bound_cond_type = "Teff";
	double F0_gauss = 1e37;
	double sigma_for_F_gauss = 5.;
	double r_gauss_cut_to_r_out = 0.01;
	double power_order = 6.;
	double kMdot_out = 2.;
	string output_dir = "data";
	bool output_fulldata = false;
	string initial_cond_shape = "power";
	string opacity_type = "Kramers";

	double Mdot_in = 0.;
	double Mdot_in_prev;
	double Mdot_out = 0.;

	{
		po::options_description desc("Allowed options");
		desc.add_options()
			( "help,h", "Produce help message" )
			( "fulldata,f", "Output files with radial structure for every computed time step. Default is output only sum.dat with integrated parameters for every time step" )
			( "alpha,a", po::value<double>(&alpha)->default_value(alpha), "Alpha parameter" )
			( "kerr,A", po::value<double>(&kerr)->default_value(kerr), "Kerr parameter of the black hole" )
			( "dilution,D", po::value<double>(&fc)->default_value(fc), "Dilution parameter"  )
			( "Mopt,m",	po::value<double>()->default_value(Mopt/GSL_CONST_CGSM_SOLAR_MASS), "Mass of optical star, solar masses" )
			( "Mx,M", po::value<double>()->default_value(Mx/GSL_CONST_CGSM_SOLAR_MASS), "Mass of central object, solar masses" )
			( "period,P", po::value<double>()->default_value(P/DAY), "Orbital period of binary system, days" )
			( "rout,R", po::value<double>()->default_value(r_out/solar_radius), "Outer radius of the disk, solar radii. If it isn't setted than it will be calculated using Mx, Mopt and period" )
			( "distance,r", po::value<double>()->default_value(Distance/kpc), "Distance to the system, kpc" )
			( "opacity,O", po::value<string>(&opacity_type)->default_value(opacity_type), "Opacity law: Kramers or OPAL" )
			( "numin,u", po::value<double>()->default_value(nu_min/keV), "Lower bound of X-ray band, keV" )
			( "numax,U", po::value<double>()->default_value(nu_max/keV), "Upper bound of X-ray band, keV" )
			( "Nx,N",	po::value<int>(&Nx)->default_value(Nx), "Size of calculation grid" )
			( "gridscale,g", po::value<string>(&grid_scale)->default_value(grid_scale), "Type of grid: log or linear" )
			( "tau,t",	po::value<double>()->default_value(tau/DAY), "Time step, days" )
			( "time,T", po::value<double>()->default_value(Time/DAY), "Computation time, days" )
			( "boundcond,B", po::value<string>(&bound_cond_type)->default_value(bound_cond_type), "Boundary movement condition, should be one of: Teff, Tirr, fourSigmaCrit, MdotOut" )
			( "Thot,H", po::value<double>(&T_min_hot_disk)->default_value(T_min_hot_disk), "Minimum photosphere temperature of the outer edge of the hot disk, degrees Kelvin. This option works only with --boundcond=Teff" )
			( "Cirr,C", po::value<double>(&C_irr)->default_value(C_irr), "Irradiation factor" )
			( "dir,d", po::value<string>(&output_dir)->default_value(output_dir), "Directory to write output files. It should exists" )
			( "F0,F", po::value<double>(&F0_gauss)->default_value(F0_gauss), "Initial viscous torque per radian on outer border of the disk, cgs" )
			( "initialcond,I", po::value<string>(&initial_cond_shape)->default_value(initial_cond_shape), "One of the available shapes of initial conditions for viscous torque F: sinusgauss, power, sinus" )
			( "powerorder,p", po::value<double>(&power_order)->default_value(power_order), "Parameter of initial condition distribution: F ~ h^poweroder. This option works only with --initialcond=power" )
		;
		po::variables_map vm;

		try {
			po::store( po::parse_command_line(ac, av, desc), vm );
			po::notify(vm);
		} catch (exception &e){
			cerr << "Error: " << e.what() << endl;
			return 1;
		}

		if ( vm.count("help") ){
			cout << desc << endl;
			return 0;
		}
		output_fulldata = vm.count("fulldata");
		Mopt = vm["Mopt"].as<double>() * GSL_CONST_CGSM_SOLAR_MASS;
		Mx = vm["Mx"].as<double>() * GSL_CONST_CGSM_SOLAR_MASS;
		P = vm["period"].as<double>() * DAY;
		Distance = vm["distance"].as<double>() * kpc;
		nu_min = vm["numin"].as<double>() * keV;
		nu_max = vm["numax"].as<double>() * keV;
		tau = vm["tau"].as<double>() * DAY;
		Time = vm["time"].as<double>() * DAY;
		if ( not vm["rout"].defaulted() ){
			r_out = vm["rout"].as<double>() * solar_radius;
		}
		else{
			r_out = r_out_func( Mx, Mopt, P );
		}
		if ( not vm["kerr"].defaulted() ){
			eta = efficiency_of_accretion(kerr);
		}
		r_in = r_in_func( Mx, kerr );
	}

	const double GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx;
	const double h_in = sqrt( GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx  * r_in );
	const double h_out = sqrt( GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx  * r_out );
	const double cosiOverD2 = cosi / Distance / Distance;

	double m_op, n_op, varkappa0, Pi_Sigma, Pi_Height, D_op, Height_exp_F, Height_exp_R, Height_coef;
	if ( opacity_type == "Kramers" ){
		// varkappa = vakappa0 * rho^1 / T^3.5
		m_op = 0.3;
		n_op = 0.8;
		varkappa0 = 5e24;

		// tau_0 = 1e3
		const double Pi1 = 6.31217;
		const double Pi2 = 0.51523;
		const double Pi3 = 1.13050;
		const double Pi4 = 0.39842;

		Pi_Sigma = pow(Pi1, 0.05) * pow(Pi2, 0.1) * pow(Pi3, 0.8) * pow(Pi4, 0.1);
		Pi_Height = pow(Pi1, 19./40.) * pow(Pi2, -0.05) * pow(Pi3, 0.1) * pow(Pi4, -0.05);

		D_op = 2.41e34 * pow(alpha, 0.8) * (Mx/GSL_CONST_CGSM_SOLAR_MASS) * pow(mu/0.5, -0.75) / Pi_Sigma * pow(varkappa0/1e22, 0.1);

		Height_exp_F = 3./20.;
		Height_exp_R = 1./8.;
		Height_coef = 0.020 * pow(1e17, -Height_exp_F) * pow(GM, -Height_exp_F/2.) * pow(1e10, -Height_exp_F/2.) * pow(alpha, -0.1) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -3./8.) * pow(mu/0.6, -3./8.) * Pi_Height * pow(varkappa0/5e24, 0.05);
	} else if ( opacity_type == "OPAL" ){
		// varkappa = vakappa0 * rho^1 / T^2.5
		m_op = 1./3.;
		n_op = 1.;
		varkappa0 = 1.5e20;

		// tau_0 = 1e3
		const double Pi1 = 5.53450;
		const double Pi2 = 0.55220;
		const double Pi3 = 1.14712;
		const double Pi4 = 0.44777;

		Pi_Sigma = pow(Pi1, 1./18.) * pow(Pi2, 1./9.) * pow(Pi3, 7./9.) * pow(Pi4, 1./9.);
		Pi_Height = pow(Pi1, 17./36.) * pow(Pi2, -1./18.) * pow(Pi3, 1./9.) * pow(Pi4, -1./18.);

		D_op = 2.12e37 * pow(alpha, 7./9.) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, 10./9.) * pow(mu/0.5, -13./18.) / Pi_Sigma * pow(varkappa0/1e22, 1./9.);
		
		Height_exp_F = 1./6.;
		Height_exp_R = 1./12.;
		Height_coef = 0.021 * pow(1e17, -Height_exp_F) * pow(GM, -Height_exp_F/2.) * pow(1e10, -Height_exp_F/2.) * pow(alpha, -1./9.) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -13./36.) * pow(mu/0.6, -13./36.) * Pi_Height * pow(varkappa0/1.5e20, 1./18.);
	} else{
		throw po::invalid_option_value(opacity_type);
	}


//	cout
//		<< "R_out = " << r_out / solar_radius << "\n"
//		<< "h_out = " << h_out << "\n"
//		<< std::flush;

	auto wunc = [m_op, n_op, D_op](
		const vector<double> &h, const vector<double> &F,
		int first, int last	
	) ->vector<double>{
		vector<double> W( first > 0 ? first : 0,  0. );
		for ( int i = first; i <= last; ++i ){
			W.push_back(
				pow(2.*M_PI*F.at(i), 1. - m_op) * pow(h.at(i), n_op) / (1. - m_op) / D_op / (2.*M_PI)
			);
		}
		return W;
	};

	// Equation from Menou et al. 1999. Sigma_cr is from their fig 8 and connected to point where Mdot is minimal. Our Sigma is 0.5 from their Sigma.
	auto Sigma_hot_disk = [alpha, Mx](double r) ->double{
		return 0.5 * 39.9 * pow(alpha/0.1, -0.80) * pow(r/1e10, 1.11) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -0.37);
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
		const double F0_sinus = 1e-6 * F0_gauss;
		const double h_cut_for_F_gauss = h_out / sqrt(r_gauss_cut_to_r_out);
		const double F_gauss_cut = F0_gauss * exp( - (h_cut_for_F_gauss-h_out)*(h_cut_for_F_gauss-h_out) / (2. * h_out*h_out/(sigma_for_F_gauss*sigma_for_F_gauss)) );
		for ( int i = 0; i < Nx; ++i ){
			double F_gauss = F0_gauss * exp( - (h.at(i)-h_out)*(h.at(i)-h_out) / (2. * h_out*h_out/(sigma_for_F_gauss*sigma_for_F_gauss)) ) - F_gauss_cut;
			F_gauss = F_gauss >= 0 ? F_gauss : 0.;
			const double F_sinus =  F0_sinus * sin( (h.at(i) - h_in) / (h_out - h_in) * M_PI / 2. );
			F.at(i) = F_gauss + F_sinus;
		}
	} else if ( initial_cond_shape == "power" ){
		for ( int i = 0; i < Nx; ++i ){
			F.at(i) = F0_gauss * pow( (h.at(i) - h_in) / (h_out - h_in), power_order );
		}
	} else if ( initial_cond_shape == "sinus" ){
		for ( int i = 0; i < Nx; ++i ){
			F.at(i) = F0_gauss * sin( (h.at(i) - h_in) / (h_out - h_in) * M_PI / 2. );
		}
	} else if ( initial_cond_shape == "sinusparabola" ){
		const double h_F0 = h_out * 0.9;
		const double delta_h = h_out - h_F0;

		F0_gauss = 1.24e13 * pow(Sigma_hot_disk(R.at(Nx-1)), 10./7.) * pow(h.at(Nx-1), 22./7.) * pow(GM, -10./7.) * pow(alpha, 8./7.);

		Mdot_out = -kMdot_out * F0_gauss / (h_F0 - h_in) * M_PI*M_PI;

		for ( int i = 0; i < Nx; ++i ){
			if ( h.at(i) < h_F0 ){
				F.at(i) = F0_gauss * sin( (h.at(i) - h_in) / (h_F0 - h_in) * M_PI / 2. );
			} else{
				F.at(i) = F0_gauss * ( 1. - kMdot_out / (h_F0-h_in) / delta_h * M_PI / 4. * (h.at(i) - h_F0)*(h.at(i) - h_F0) );
			}
		}
	} else{
		throw po::invalid_option_value(initial_cond_shape);
	}	

	ofstream output_sum( output_dir + "/sum.dat" );
	output_sum << "#t	Mdot	Lx	H2R	Rhot2Rout	Tphout	kxout Qiir2Qvisout	Qirr2Qvisout_analyt	mB	mV" << "\n";
	output_sum << "# r_out = " << r_out << "\n";
	output_sum << "#";
	for ( int i = 0; i < ac; ++i ){
		output_sum << " " << av[i];
	}
	output_sum << endl;

	for( double t = 0.; t <= Time; t += tau ){
		// cout << t/DAY << endl;

		vector<double> W(Nx, 0.), Tph(Nx, 0.), Tph_vis(Nx, 0.), Tph_X(Nx, 0.), Tirr(Nx,0.), Sigma(Nx, 0.), Height(Nx, 0.);
		
		try{
			nonlenear_diffusion_nonuniform_1_2 (tau, eps, 0., Mdot_out/(2.*M_PI), wunc, h, F);
			W = wunc(h, F, 1, Nx-1);
		} catch (runtime_error er){
			cout << er.what() << endl;
			break;
		}

		Mdot_in_prev = Mdot_in;
		Mdot_in = 2.*M_PI * ( F.at(1) - F.at(0) ) / ( h.at(1) - h.at(0) );

		for ( int i = 1; i < Nx; ++i ){
			Sigma.at(i) = 0.5 * W.at(i) * GM*GM / (2. * pow(h.at(i), 3.));
			Height.at(i) = R.at(i) * Height_coef * pow(2.*M_PI*F.at(i), Height_exp_F) * pow(R.at(i)/1e10, Height_exp_R - Height_exp_F/2.);
			Tph_vis.at(i) = GM * pow(h.at(i), -1.75) * pow( 0.75 * F.at(i) / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph_X.at(i) = fc * T_GR( R.at(i), 0., Mx, Mdot_in, R.front() );

			// k_x = k_irr * pow(Height.at(ii) / R.at(ii), 2.);
			const double Qx = C_irr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (4.*M_PI * R.at(i)*R.at(i));
			Tirr.at(i) = pow( Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph.at(i) = pow( pow(Tph_vis.at(i), 4.) + Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
		}
		
		const double Lx = Luminosity( R, Tph_X, nu_min, nu_max, 100 ) / pow(fc, 4.);

		double k_x = 0.;
		int ii = Nx;
		if (bound_cond_type == "MdotOut"){
			Mdot_out = - kMdot_out * Mdot_in;
			do{
				ii--;
			} while( Sigma.at(ii) < Sigma_hot_disk(R[ii]) );

		} else if (bound_cond_type == "fourSigmaCrit"){
			do{
				ii--;
				// Equation from Menou et al. 1999. Factor 4 is from their fig 8 and connected to point where Mdot = 0. Our Sigma is 0.5 from their Sigma.
			} while( Sigma.at(ii) <  4 * Sigma_hot_disk(R[ii]) );
		} else if ( bound_cond_type == "Teff" ){
			do{
				ii--;
			} while( Tph.at(ii) < T_min_hot_disk );
		} else if (  bound_cond_type == "Tirr" ){
			if ( Mdot_in >= Mdot_in_prev ){
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
		
		const double mB = -2.5 * log10( I_lambda(R, Tph, lambdaB) * cosiOverD2 / irr0B );
		const double mV = -2.5 * log10( I_lambda(R, Tph, lambdaV) * cosiOverD2 / irr0V );

		if (output_fulldata){
			ostringstream filename;
			filename << output_dir << "/" << static_cast<int>(t/tau) << ".dat";
			ofstream output( filename.str() );
			output << "#h   F   Sigma   W   R   Tph_vis Height	Tph" << "\n";
			output << "# Time = " << t / DAY << " Mdot_in = " << Mdot_in << endl;
			for ( int i = 1; i < Nx; ++i ){
				output		<< h.at(i)
					<< "\t" << F.at(i)
					<< "\t" << Sigma.at(i)
					<< "\t" << W.at(i)
					<< "\t" << R.at(i)
					<< "\t" << Tph_vis.at(i)
					<< "\t" << Height.at(i)
					<< "\t" << Tph.at(i)
					<< endl;
			}
		}

		output_sum		<< t / DAY
				<< "\t" << Mdot_in
				<< "\t" << Lx
				<< "\t" << Height.at(Nx-1) / R.at(Nx-1)
				<< "\t" << R.at(Nx-1) / r_out
				<< "\t" << Tph.at(Nx-1)
				<< "\t" << C_irr
				<< "\t" << pow( Tirr.at(Nx-1) / Tph_vis.at(Nx-1), 4. )
				<< "\t" << 4./3. * eta * C_irr * R.at(Nx-1) / (2. * GM / GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_SPEED_OF_LIGHT)
				<< "\t" << mB
				<< "\t" << mV
				<< endl;

		if ( ii < Nx-1 ){
			Nx = ii+1;
			// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
			h.resize(Nx);
		}
	}

	return 0;
}
