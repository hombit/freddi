#include "freddi.hpp"

//#include <algorithm>
//#include <boost/any.hpp>
#include <cmath>
#include <exception>
#include <fstream>
//#include <iostream>
//#include <limits>
//#include <sstream>
//#include <string>
#include <vector>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "constants.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"
#include "spectrum.hpp"


//namespace po = boost::program_options;
//namespace odeint = boost::numeric::odeint;
using namespace std;
//using namespace std::placeholders;


void freddi(const FreddiArguments& args) {
	double Mdot_in = 0.;
	double Mdot_in_prev;
	double Mdot_out = 0.;
	double Mdisk = 0.;

	const double GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args.basic->Mx;
	const double eta = disk_orbit::efficiency_of_accretion(args.basic->kerr);
	const double h_in = args.basic->h(args.basic->rin);
	const double h_out = args.basic->h(args.basic->rout);
	const double cosi = std::cos(args.flux->inclination / 180 * M_PI);
	const double cosiOverD2 = cosi / (args.flux->distance * args.flux->distance);
	const OpacityRelated* oprel = args.disk->oprel.get();
	auto Nx = args.calc->Nx;

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
	auto Sigma_hot_disk = [args](double r) ->double{
		return 39.9 * pow(args.basic->alpha/0.1, -0.80) * pow(r/1e10, 1.11) * pow(args.basic->Mx/GSL_CONST_CGSM_SOLAR_MASS, -0.37);
	};

	vector<double> h(Nx), R(Nx);
	for ( int i = 0; i < Nx; ++i ){
		if (args.calc->gridscale == "log") {
			h[i] = h_in * pow( h_out/h_in, i/(Nx-1.) );
		} else if (args.calc->gridscale == "linear") {
			h[i] = h_in + (h_out - h_in) * i/(Nx-1.);
		} else{
			throw logic_error("Wrong gridscale");
		}
		R[i] = h[i] * h[i] / GM;
	}

	vector<double> F(Nx);
	if (args.disk->initialcond == "power" or args.disk->initialcond == "powerF") {
		for ( int i = 0; i < Nx; ++i ){
			F[i] = args.disk->F0 * pow((h[i] - h_in) / (h_out - h_in), args.disk->powerorder);
		}
	} else if (args.disk->initialcond == "powerSigma") {
		for ( int i = 0; i < Nx; ++i ){
			const double Sigma_to_Sigmaout = pow( (h[i] - h_in) / (h_out - h_in), args.disk->powerorder );
			F[i] = args.disk->F0 * pow( h[i] / h_out, (3. - oprel->n) / (1. - oprel->m) ) * pow( Sigma_to_Sigmaout, 1. / (1. - oprel->m) );
		}
	} else if(args.disk->initialcond == "quasistat") {
		for ( int i = 0; i < Nx; ++i ){
			const double xi_LS2000 = h[i] / h_out;
			F[i] = args.disk->F0 * oprel->f_F(xi_LS2000) * (1. - h_in / h[i]) / (1. - h_in / h_out);
		}
	} else if(args.disk->initialcond == "gaussF"){
		for ( int i = 0; i < Nx; ++i ){
			const double xi = (h.at(i) - h_in) / (h_out - h_in);
			F[i] = args.disk->F0 * exp( -(xi - args.disk->gaussmu)*(xi - args.disk->gaussmu) / (2. * args.disk->gausssigma*args.disk->gausssigma) );
		}
	} else{
		throw logic_error("Wrong initialcond");
	}

	const std::string output_filepath(args.general->dir + "/" + args.general->prefix + ".dat");
	ofstream output_sum(output_filepath, std::ofstream::out);
	output_sum << "#t    Mdot Mdisk Rhot Cirrout H2R   Teffout Tirrout Qiir2Qvisout Lx    mU  mB  mV  mR  mI  mJ ";
	for ( int i = 0; i < args.flux->lambdas.size(); ++i ){
		output_sum << " Fnu" << i;
		for ( double j = 0; j < 9 - log10(i+0.1); ++j ){
			output_sum << " ";
		}
	}
	output_sum << "\n";
	output_sum << "#days g/s  g     Rsun float   float K       K       float        erg/s mag mag mag mag mag mag";
	for ( int i = 0; i < args.flux->lambdas.size(); ++i ){
		output_sum << " erg/s/cm^2/Hz";
	}
	output_sum << "\n";
	// TODO: move file output somewhere
//	for ( const auto &it : po_vm ){
//		auto &value = it.second.value();
//		if ( auto v = boost::any_cast<uint32_t>(&value) ){
//			output_sum << "# "
//			           << it.first.c_str()
//			           << "="
//			           << *v
//					   << "\n";
//		} else if ( auto v = boost::any_cast<string>(&value) ){
//			output_sum << "# "
//			           << it.first.c_str()
//			           << "="
//			           << *v
//					   << "\n";
//		} else if ( auto v = boost::any_cast<double>(&value) ){
//			output_sum << "# "
//			           << it.first.c_str()
//			           << "="
//			           << *v
//					   << "\n";
//		} else if ( auto v = boost::any_cast<unsigned int>(&value) ){
//			output_sum << "# "
//			           << it.first.c_str()
//			           << "="
//			           << *v
//					   << "\n";
//		} else if ( auto v = boost::any_cast< vector<double> >(&value) ){
//			for ( int i = 0; i < v->size(); ++i ){
//				output_sum << "# "
//						   << it.first.c_str()
//						   << "="
//						   << v->at(i)
//						   << "  # "
//						   << i
//						   << "\n";
//				}
//		} else {
//			output_sum << "error\n";
//			// throw po::invalid_option_value(it.first.c_str());
//		}
//	}
//	if ( po_vm.count("rout") == 0 ){
//		output_sum << "# --rout hadn't been specified, tidal radius " << r_out/solar_radius << " Rsun was used" << std::endl;
//	}
	output_sum << flush;

	vector<double> W(Nx, 0.);
	W = wunc(h, F, 1, Nx-1);

	for( int i_t = 0; i_t <= args.calc->time/args.calc->tau + 0.001; ++i_t ){
		const double t = i_t * args.calc->tau;
		// cout << t/DAY << endl;

		vector<double> Tph(Nx, 0.), Tph_vis(Nx, 0.), Tph_X(Nx, 0.), Tirr(Nx,0.), Sigma(Nx, 0.), Height(Nx, 0.);

		Mdot_in_prev = Mdot_in;
		Mdot_in = ( F.at(1) - F.at(0) ) / ( h.at(1) - h.at(0) );

		double C_irr;
		for ( int i = 1; i < Nx; ++i ){
			Sigma[i] = W[i] * GM*GM / ( 4.*M_PI *  pow(h[i], 3.) );
			Height[i] = oprel->Height(R[i], F[i]);
			Tph_vis[i] = GM * pow(h[i], -1.75) * pow( 3. / (8.*M_PI) * F[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph_X[i] = args.flux->colourfactor * T_GR( R[i], args.basic->kerr, args.basic->Mx, Mdot_in, R.front() );

			double Qx;
			if (args.irr->irrfactortype == "const") {
				C_irr = args.irr->Cirr;
				Qx = args.irr->Cirr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (4.*M_PI * R[i]*R[i]);
			} else if (args.irr->irrfactortype == "square") {
				C_irr = args.irr->Cirr * (Height[i] / R[i]) * (Height[i] / R[i]);
				Qx = C_irr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (4.*M_PI * R[i]*R[i]);
			} else{
				throw logic_error("Wrong irrfactor");
			}
			Tirr[i] = pow( Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
			Tph[i] = pow( pow(Tph_vis[i], 4.) + Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25 );
		}

		const double Lx = Luminosity(R, Tph_X, args.flux->emin, args.flux->emax, 100) / pow(args.flux->colourfactor, 4.);

		const double mU = -2.5 * log10( I_lambda(R, Tph, lambdaU) * cosiOverD2 / irr0U );
		const double mB = -2.5 * log10( I_lambda(R, Tph, lambdaB) * cosiOverD2 / irr0B );
		const double mV = -2.5 * log10( I_lambda(R, Tph, lambdaV) * cosiOverD2 / irr0V );
		const double mR = -2.5 * log10( I_lambda(R, Tph, lambdaR) * cosiOverD2 / irr0R );
		const double mI = -2.5 * log10( I_lambda(R, Tph, lambdaI) * cosiOverD2 / irr0I );
		const double mJ = -2.5 * log10( I_lambda(R, Tph, lambdaJ) * cosiOverD2 / irr0J );

		Mdisk = 0.;
		for ( int i = 0; i < Nx; ++i ){
			double stepR;
			if ( i == 0              ) stepR = R[i+1] - R[i  ];
			if ( i == Nx-1           ) stepR = R[i  ] - R[i-1];
			if ( i > 1 and i < Nx-1  ) stepR = R[i+1] - R[i-1];
			Mdisk += 0.5 * Sigma[i] * 2.*M_PI * R[i] * stepR;
		}

		if (args.general->fulldata){
			ostringstream filename;
			filename << args.general->dir << "/" << args.general->prefix << "_" << i_t << ".dat";
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
		for ( auto &lambda : args.flux->lambdas ){
			output_sum
			    << "\t" << I_lambda(R, Tph, lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * cosiOverD2;
		}
		output_sum      << endl;


		try{
			nonlenear_diffusion_nonuniform_1_2 (args.calc->tau, args.calc->eps, 0., Mdot_out, wunc, h, F);
			W = wunc(h, F, 1, Nx-1);
		} catch (runtime_error er){
			cout << er.what() << endl;
			break;
		}


		int ii = Nx;
//		if (bound_cond_type == "MdotOut"){
//			Mdot_out = - kMdot_out * Mdot_in;
//			do{
//				ii--;
//			} while( Sigma.at(ii) < Sigma_hot_disk(R[ii]) );
//
//		} else if (bound_cond_type == "fourSigmaCrit"){
//			do{
//				ii--;
//				 Equation from Menou et al. 1999. Factor 4 is from their fig 8 and connected to point where Mdot = 0.
//			} while( Sigma.at(ii) <  4 * Sigma_hot_disk(R[ii]) );
//		} else
		if (args.disk->boundcond == "Teff") {
			do{
				ii--;
			} while( Tph.at(ii) < args.disk->Thot );
		} else if (args.disk->boundcond == "Tirr") {
			if ( Mdot_in >= Mdot_in_prev && (args.disk->initialcond == "power" || args.disk->initialcond == "sinusgauss") ){
				do{
					ii--;
				} while( Tph.at(ii) < args.disk->Thot );
			} else{
				do{
					ii--;
				} while( Tirr.at(ii) < args.disk->Thot );
			}
		} else{
			throw logic_error("Wrong boundcond");
		}

		if ( ii < Nx-1 ){
			Nx = ii+1;
			// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
			h.resize(Nx);
		}
	}
}
