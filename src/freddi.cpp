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
		cosi(std::cos(args_.flux->inclination / 180 * M_PI)),
		cosiOverD2(std::cos(args_.flux->inclination / 180 * M_PI) / (args_.flux->distance * args_.flux->distance)),
		oprel(args_.disk->oprel.get()),
		wunc(std::bind(&Freddi::wunction, this, _1, _2, _3, _4)),
		state_(new FreddiState(initializeState())) {
	calculateRadialStructure();
}

FreddiState Freddi::initializeState() {
	FreddiState state(this);

	const double h_in = args->basic->h(args->basic->rin);
	const double h_out = args->basic->h(args->basic->rout);

	for (int i = 0; i < state.Nx; ++i) {
		if (args->calc->gridscale == "log") {
			state.h[i] = h_in * pow(h_out / h_in, i / (state.Nx - 1.));
		} else if (args->calc->gridscale == "linear") {
			state.h[i] = h_in + (h_out - h_in) * i / (state.Nx - 1.);
		} else {
			throw std::logic_error("Wrong gridscale");
		}
		state.R[i] = state.h[i] * state.h[i] / GM;
	}

	if (args->disk->initialcond == "power" or args->disk->initialcond == "powerF") {
		for (int i = 0; i < state.Nx; ++i) {
			state.F[i] = args->disk->F0 * pow((state.h[i] - h_in) / (h_out - h_in), args->disk->powerorder);
		}
	} else if (args->disk->initialcond == "powerSigma") {
		for (int i = 0; i < state.Nx; ++i) {
			const double Sigma_to_Sigmaout = pow((state.h[i] - h_in) / (h_out - h_in), args->disk->powerorder);
			state.F[i] = args->disk->F0 * pow(state.h[i] / h_out, (3. - oprel->n) / (1. - oprel->m)) *
				   pow(Sigma_to_Sigmaout, 1. / (1. - oprel->m));
		}
	} else if (args->disk->initialcond == "quasistat") {
		for (int i = 0; i < state.Nx; ++i) {
			const double xi_LS2000 = state.h[i] / h_out;
			state.F[i] = args->disk->F0 * oprel->f_F(xi_LS2000) * (1. - h_in / state.h[i]) / (1. - h_in / h_out);
		}
	} else if (args->disk->initialcond == "gaussF") {
		for (int i = 0; i < state.Nx; ++i) {
			const double xi = (state.h[i] - h_in) / (h_out - h_in);
			state.F[i] = args->disk->F0 * exp(-(xi - args->disk->gaussmu) * (xi - args->disk->gaussmu) /
										(2. * args->disk->gausssigma * args->disk->gausssigma));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}

	return state;
}

void Freddi::calculateRadialStructure() {
	state_->W = wunc(state_->h, state_->F, 1, state_->Nx - 1);

	Mdot_in_prev = state_->Mdot_in;
	state_->Mdot_in = (state_->F.at(1) - state_->F.at(0)) / (state_->h.at(1) - state_->h.at(0));

	for (int i = 1; i < state_->Nx; ++i) {
		state_->Sigma[i] = state_->W[i] * GM * GM / (4. * M_PI * pow(state_->h[i], 3.));
		state_->Height[i] = oprel->Height(state_->R[i], state_->F[i]);
		state_->Tph_vis[i] = GM * std::pow(state_->h[i], -1.75) *
				std::pow(3. / (8. * M_PI) * state_->F[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		state_->Tph_X[i] = args->flux->colourfactor * T_GR(
				state_->R[i],
				args->basic->kerr,
				args->basic->Mx,
				state_->Mdot_in,
				state_->R.front());

		double Qx;
		if (args->irr->irrfactortype == "const") {
			state_->Cirr[i] = args->irr->Cirr;
			Qx = state_->Cirr[i] * eta * state_->Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT /
				 (4. * M_PI * state_->R[i] * state_->R[i]);
		} else if (args->irr->irrfactortype == "square") {
			state_->Cirr[i] = args->irr->Cirr * (state_->Height[i] / state_->R[i]) * (state_->Height[i] / state_->R[i]);
			Qx = state_->Cirr[i] * eta * state_->Mdot_in * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT /
				 (4. * M_PI * state_->R[i] * state_->R[i]);
		} else {
			throw std::logic_error("Wrong irrfactor");
		}
		state_->Tirr[i] = pow(Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		state_->Tph[i] = pow(pow(state_->Tph_vis[i], 4.) + Qx / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
	}

	state_->Lx = Luminosity(state_->R, state_->Tph_X, args->flux->emin, args->flux->emax, 100) / pow(args->flux->colourfactor, 4.);
}

void Freddi::next() {
	state_.reset(new FreddiState(*state_));

	nonlenear_diffusion_nonuniform_1_2(args->calc->tau, args->calc->eps, 0., state_->Mdot_out, wunc, state_->h, state_->F);

	state_->increase_t(args->calc->tau);
	calculateRadialStructure();
	truncateOuterRadius();
}


void Freddi::truncateOuterRadius() {
	int ii = state_->Nx;
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
		} while( state_->Tph.at(ii) < args->disk->Thot );
	} else if (args->disk->boundcond == "Tirr") {
		if ( state_->Mdot_in >= Mdot_in_prev && (args->disk->initialcond == "power" || args->disk->initialcond == "sinusgauss") ){
			do{
				ii--;
			} while( state_->Tph.at(ii) < args->disk->Thot );
		} else{
			do{
				ii--;
			} while( state_->Tirr.at(ii) < args->disk->Thot );
		}
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii < state_->Nx-1 ){
		state_->Nx = ii+1;
		// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
		state_->h.resize(state_->Nx);
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



FreddiState::FreddiState(const Freddi *freddi):
		freddi(freddi),
		Nx(freddi->args->calc->Nx),
		h(Nx),
		R(Nx),
		F(Nx),
		W(Nx, 0.),
		Tph(Nx, 0.),
		Tph_vis(Nx, 0.),
		Tph_X(Nx, 0.),
		Tirr(Nx, 0.),
		Cirr(Nx, 0.),
		Sigma(Nx, 0.),
		Height(Nx, 0.) {}


double FreddiState::integrate(const vecd& values) const {
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


void FreddiState::increase_t(double tau) {
	t += tau;
	i_t += 1;
}
