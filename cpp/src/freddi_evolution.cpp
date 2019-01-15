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


FreddiEvolution::FreddiEvolution(const FreddiArguments &args):
		state_(FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4))) {}


void FreddiEvolution::step(const double tau) {
	state_.before_step(tau);
	truncateInnerRadius();
	nonlinear_diffusion_nonuniform_wind_1_2(
			state_.args.calc->tau, state_.args.calc->eps,
			state_.F_in(), state_.Mdot_out(),
			state_.wind_->A(), state_.wind_->B(), state_.wind_->C(),
			state_.wunc,
			state_.h(), state_.F_);
	truncateOuterRadius();
}


void FreddiEvolution::truncateOuterRadius() {
	if (state_.Mdot_in() > state_.Mdot_in_prev()) {
		return;
	}

	auto ii = state_.Nx();
	if (state_.args.disk->boundcond == "Teff") {
		do {
			ii--;
		} while( state_.Tph().at(ii) < state_.args.disk->Thot );
	} else if (state_.args.disk->boundcond == "Tirr") {
		do {
			ii--;
		} while( state_.Tirr().at(ii) < state_.args.disk->Thot );
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii < state_.Nx_ - 1 ){
		state_.Nx_ = ii+1;
		// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
		state_.h_.resize(state_.Nx_);
		state_.R_.resize(state_.Nx_);
		state_.F_.resize(state_.Nx_);
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t first, size_t last) const {
	vecd W(last + 1, 0.);
	for ( size_t i = first; i <= last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - state_.oprel.m) * pow(h[i], state_.oprel.n) / (1. - state_.oprel.m) / state_.oprel.D;
	}
	return W;
};



FreddiNeutronStarEvolution::FreddiNeutronStarEvolution(const FreddiNeutronStarArguments &args):
		FreddiEvolution(args),
		args_ns(args.ns.get()),
		xi_pow_minus_7_2(std::pow(xi, -3.5)),
		R_m_min(std::max(args.ns->Rx, args.basic->rin)),
		mu_magn(0.5 * args.ns->Bx * args.ns->Rx*args.ns->Rx*args.ns->Rx),
		R_dead(std::cbrt(mu_magn*mu_magn / args.ns->Fdead)),
		R_cor(std::cbrt(state_.GM / (4 * M_PI*M_PI * args.ns->freqx*args.ns->freqx))) {}


void FreddiNeutronStarEvolution::truncateInnerRadius() {
	if (args_ns->Fdead <= 0.) {
		return;
	}
	if ( state_.Mdot_in() > state_.Mdot_in_prev() ) {
		return;
	}

	const double R_alfven = args_ns->epsilonAlfven *
			std::pow(mu_magn*mu_magn*mu_magn*mu_magn / (state_.Mdot_in()*state_.Mdot_in() * state_.GM), 1./7.);
	double R_m = std::max(R_m_min, R_alfven);
	R_m = std::min(R_m, R_dead);
	size_t ii;
	for (ii = 0; ii < state_.Nx() - 1; ii++) {
		if (state_.R().at(ii+1) > R_m){
			break;
		}
	}
	if (ii >= state_.Nx() - 2) {
		throw std::runtime_error("R_in > R_out");
	}
	R_m = state_.R().at(ii);

	double new_F_in;
	if (R_m < R_cor) {
		const double n_ws = 1 - k_t * xi_pow_minus_7_2 * std::pow(R_m / R_cor, 3);
		new_F_in = (1 - n_ws) * state_.Mdot_in() * std::sqrt(state_.GM * R_m);
	} else {
		new_F_in = args_ns->Fdead * std::pow(R_dead / R_m, 3);
	}
	state_.F_in_ = state_.F_[0] = new_F_in;

	if (ii > 0) {
		state_.Nx_ -= ii;
		state_.h_.erase(state_.h_.begin(), state_.h_.begin() + ii);
		state_.R_.erase(state_.R_.begin(), state_.R_.begin() + ii);
		state_.F_.erase(state_.F_.begin(), state_.F_.begin() + ii);
	}
}
