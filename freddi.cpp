#include "freddi.hpp"

#include <cmath>
#include <exception>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "constants.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"

using namespace std::placeholders;

Freddi::Freddi(const FreddiArguments &args_):
		args(&args_),
		GM(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args_.basic->Mx),
		eta(disk_orbit::efficiency_of_accretion(args_.basic->kerr)),
		h_in(args_.basic->h(args_.basic->rin)),
		h_out(args_.basic->h(args_.basic->rout)),
		cosi(std::cos(args_.flux->inclination / 180 * M_PI)),
		cosiOverD2(std::cos(args_.flux->inclination / 180 * M_PI) / (args_.flux->distance * args_.flux->distance)),
		oprel(args_.disk->oprel.get()),
		wunc(std::bind(&Freddi::wunction, this, _1, _2, _3, _4)),
		Nx(args_.calc->Nx),
		h(args_.calc->Nx),
		R(args_.calc->Nx),
		F(args_.calc->Nx),
		W(args_.calc->Nx, 0.),
		Tph(args_.calc->Nx, 0.),
		Tph_vis(args_.calc->Nx, 0.),
		Tph_X(args_.calc->Nx, 0.),
		Tirr(args_.calc->Nx, 0.),
		Sigma(args_.calc->Nx, 0.),
		Height(args_.calc->Nx, 0.) {
	initializeRadialStructure();
}

void Freddi::initializeRadialStructure() {
	for (int i = 0; i < Nx; ++i) {
		if (args->calc->gridscale == "log") {
			h[i] = h_in * pow(h_out / h_in, i / (Nx - 1.));
		} else if (args->calc->gridscale == "linear") {
			h[i] = h_in + (h_out - h_in) * i / (Nx - 1.);
		} else {
			throw std::logic_error("Wrong gridscale");
		}
		R[i] = h[i] * h[i] / GM;
	}

	if (args->disk->initialcond == "power" or args->disk->initialcond == "powerF") {
		for (int i = 0; i < Nx; ++i) {
			F[i] = args->disk->F0 * pow((h[i] - h_in) / (h_out - h_in), args->disk->powerorder);
		}
	} else if (args->disk->initialcond == "powerSigma") {
		for (int i = 0; i < Nx; ++i) {
			const double Sigma_to_Sigmaout = pow((h[i] - h_in) / (h_out - h_in), args->disk->powerorder);
			F[i] = args->disk->F0 * pow(h[i] / h_out, (3. - oprel->n) / (1. - oprel->m)) *
				   pow(Sigma_to_Sigmaout, 1. / (1. - oprel->m));
		}
	} else if (args->disk->initialcond == "quasistat") {
		for (int i = 0; i < Nx; ++i) {
			const double xi_LS2000 = h[i] / h_out;
			F[i] = args->disk->F0 * oprel->f_F(xi_LS2000) * (1. - h_in / h[i]) / (1. - h_in / h_out);
		}
	} else if (args->disk->initialcond == "gaussF") {
		for (int i = 0; i < Nx; ++i) {
			const double xi = (h.at(i) - h_in) / (h_out - h_in);
			F[i] = args->disk->F0 * exp(-(xi - args->disk->gaussmu) * (xi - args->disk->gaussmu) /
										(2. * args->disk->gausssigma * args->disk->gausssigma));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}

	calculateRadialStructure();
}

void Freddi::calculateRadialStructure() {
	W = wunc(h, F, 1, Nx - 1);

	Mdot_in_prev = Mdot_in;
	Mdot_in = (F.at(1) - F.at(0)) / (h.at(1) - h.at(0));

	for (int i = 1; i < Nx; ++i) {
		Sigma[i] = W[i] * GM * GM / (4. * M_PI * pow(h[i], 3.));
		Height[i] = oprel->Height(R[i], F[i]);
		Tph_vis[i] =
				GM * pow(h[i], -1.75) * pow(3. / (8. * M_PI) * F[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		Tph_X[i] = args->flux->colourfactor * T_GR(R[i], args->basic->kerr, args->basic->Mx, Mdot_in, R.front());

		double Qx;
		if (args->irr->irrfactortype == "const") {
			const double C_irr = args->irr->Cirr;
			Qx = args->irr->Cirr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT /
				 (4. * M_PI * R[i] * R[i]);
		} else if (args->irr->irrfactortype == "square") {
			const double C_irr = args->irr->Cirr * (Height[i] / R[i]) * (Height[i] / R[i]);
			Qx = C_irr * eta * Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT /
				 (4. * M_PI * R[i] * R[i]);
		} else {
			throw std::logic_error("Wrong irrfactor");
		}
		Tirr[i] = pow(Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		Tph[i] = pow(pow(Tph_vis[i], 4.) + Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
	}

	Lx = Luminosity(R, Tph_X, args->flux->emin, args->flux->emax, 100) / pow(args->flux->colourfactor, 4.);
}

void Freddi::next() {
	t += args->calc->tau;

	nonlenear_diffusion_nonuniform_1_2(args->calc->tau, args->calc->eps, 0., Mdot_out, wunc, h, F);
	W = wunc(h, F, 1, Nx - 1);

	calculateRadialStructure();
	truncateOuterRadius();
}


void Freddi::truncateOuterRadius() {
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
	if (args->disk->boundcond == "Teff") {
		do{
			ii--;
		} while( Tph.at(ii) < args->disk->Thot );
	} else if (args->disk->boundcond == "Tirr") {
		if ( Mdot_in >= Mdot_in_prev && (args->disk->initialcond == "power" || args->disk->initialcond == "sinusgauss") ){
			do{
				ii--;
			} while( Tph.at(ii) < args->disk->Thot );
		} else{
			do{
				ii--;
			} while( Tirr.at(ii) < args->disk->Thot );
		}
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii < Nx-1 ){
		Nx = ii+1;
		// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
		h.resize(Nx);
	}
}

vecd Freddi::wunction(const vecd &h, const vecd &F, int first, int last) const {
	vecd W(first > 0 ? first : 0,  0.);
	for ( int i = first; i <= last; ++i ){
		W.push_back(pow(F[i], 1. - oprel->m) * pow(h[i], oprel->n) / (1. - oprel->m) / oprel->D);
	}
	return W;
};


// Equation from Lasota, Dubus, Kruk A&A 2008, Menou et al. 1999. Sigma_cr is from their fig 8 and connected to point where Mdot is minimal.
double Freddi::Sigma_hot_disk(double r) const {
	return 39.9 * pow(args->basic->alpha/0.1, -0.80) * pow(r/1e10, 1.11) * pow(args->basic->Mx/GSL_CONST_CGSM_SOLAR_MASS, -0.37);
};


double Freddi::integrate(const vecd& values) const {
	double integral = 0.;
	size_t N = std::min(values.size(), R.size());
	double stepR;
	for ( int i = 0; i < N; ++i ){
		if ( i == 0            ) stepR = R[i+1] - R[i  ];
		if ( i == N-1          ) stepR = R[i  ] - R[i-1];
		if ( i > 1 and i < N-1 ) stepR = R[i+1] - R[i-1];
		integral += 0.5 * values[i] * 2.*M_PI * R[i] * stepR;
	}
	return integral;
}

