#include "freddi_evolution.hpp"

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
		FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)) {}


void FreddiEvolution::step(const double tau) {
	FreddiState::step(tau);
	truncateInnerRadius();
	nonlinear_diffusion_nonuniform_wind_1_2(
			args.calc->tau, args.calc->eps,
			F_in(), Mdot_out(),
			windA(), windB(), windC(),
			wunc,
			h(), F_,
			first(), last());
	truncateOuterRadius();
}


void FreddiEvolution::truncateOuterRadius() {
	if (Mdot_in() > Mdot_in_prev()) {
		return;
	}

	auto ii = last() + 1;
	if (args.disk->boundcond == "Teff") {
		do {
			ii--;
		} while( Tph().at(ii) < args.disk->Thot );
	} else if (args.disk->boundcond == "Tirr") {
		do {
			ii--;
		} while( Tirr().at(ii) < args.disk->Thot );
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii <= last() - 1 ){
		last_ = ii;
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t _first, size_t _last) const {
	vecd W(_last + 1, 0.);
	for ( size_t i = _first; i <= _last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - oprel.m) * pow(h[i], oprel.n) / (1. - oprel.m) / oprel.D;
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
		R_cor(std::cbrt(GM / (4 * M_PI*M_PI * args.ns->freqx*args.ns->freqx))),
		inverse_beta(args.ns->inversebeta),
		Fmagn(initialize_Fmagn()),
		dFmagn_dh(initialize_dFmagn_dh()),
		d2Fmagn_dh2(initialize_d2Fmagn_dh2()) {}


const vecd FreddiNeutronStarEvolution::initialize_Fmagn() {
	vecd Fmagn_(Nx);
	const double k = inverse_beta * mu_magn*mu_magn / 3.;
	double brackets;
	for (size_t i = 0; i < Nx; i++) {
		if (R_[i] < R_cor) {
			brackets = -1. + 2. * std::pow(R_[i] / R_cor, 1.5) - 2./3. * std::pow(R_[i] / R_cor, 3);
		} else {
			brackets = 1. - 2./3. * std::pow(R_cor / R_[i], 1.5);
		}
		Fmagn_[i] = k / (R_[i]*R_[i]*R_[i]) * brackets;
	}
	return Fmagn_;
}


const vecd FreddiNeutronStarEvolution::initialize_dFmagn_dh() {
	vecd dFmagn_dh_(Nx);
	const double k = inverse_beta * 2*mu_magn*mu_magn * GM*GM*GM;
	double brackets;
	for (size_t i = 0; i < Nx; i++) {
		if (R_[i] < R_cor) {
			brackets = 1. - std::pow(R_[i] / R_cor, 1.5);
		} else {
			brackets = -1. + std::pow(R_cor / R_[i], 1.5);
		}
		dFmagn_dh_[i] = k / std::pow(h_[i], 7) * brackets;
	}
	return dFmagn_dh_;
}


const vecd FreddiNeutronStarEvolution::initialize_d2Fmagn_dh2() {
	vecd d2Fmagn_dh2_(Nx);
	const double k = inverse_beta * 2*mu_magn*mu_magn / GM;
	double brackets;
	for (size_t i = 0; i < Nx; i++) {
		if (R_[i] < R_cor) {
			brackets = -7. + 4. * std::pow(R_[i] / R_cor, 1.5);
		} else {
			brackets = 7. - 10. * std::pow(R_cor / R_[i], 1.5);
		}
		d2Fmagn_dh2_[i] = k / std::pow(R_[i], 4) * brackets;
	}
	return d2Fmagn_dh2_;
}


void FreddiNeutronStarEvolution::truncateInnerRadius() {
	if (args_ns->Fdead <= 0.) {
		return;
	}
	if ( Mdot_in() > Mdot_in_prev() ) {
		return;
	}

	const double R_alfven = args_ns->epsilonAlfven *
			std::pow(mu_magn*mu_magn*mu_magn*mu_magn / (Mdot_in()*Mdot_in() * GM), 1./7.);
	double R_m = std::max(R_m_min, R_alfven);
	R_m = std::min(R_m, R_dead);
	size_t ii;
	for (ii = first_; ii <= last_ - 2; ii++) {
		if (R_[ii + 1] > R_m){
			break;
		}
	}
	if (ii >= last_ - 2) {
		throw std::runtime_error("R_in > R_out");
	}
	R_m = R_[ii];
	first_ = ii;

	double new_F_in = 0;
	if (inverse_beta <= 0.) {
		if (R_m < R_cor) {
			const double n_ws = 1 - k_t * xi_pow_minus_7_2 * std::pow(R_m / R_cor, 3);
			new_F_in = (1 - n_ws) * Mdot_in() * std::sqrt(GM * R_m);
		} else {
			new_F_in = args_ns->Fdead * std::pow(R_dead / R_m, 3);
		}
	}
	F_in_ = F_[first_] = new_F_in;
}


double FreddiNeutronStarEvolution::Mdot_in() const {
	const double dF_dh = (F_[first_ + 1] - F_[first_]) / (h_[first_ + 1] - h_[first_]);
	return dF_dh + dFmagn_dh[first_];
}


vecd FreddiNeutronStarEvolution::windC() const {
	auto C = FreddiEvolution::windC();
	for (size_t i = 0; i < C.size(); i++) {
		C[i] += d2Fmagn_dh2[i];
	}
	return C;
}
