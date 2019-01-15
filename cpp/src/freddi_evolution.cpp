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
			h(), F_);
	truncateOuterRadius();
}


void FreddiEvolution::truncateOuterRadius() {
	if (Mdot_in() > Mdot_in_prev()) {
		return;
	}

	auto ii = Nx();
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

	if ( ii < Nx_ - 1 ){
		Nx_ = ii+1;
		// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
		h_.resize(Nx_);
		R_.resize(Nx_);
		F_.resize(Nx_);
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t first, size_t last) const {
	vecd W(last + 1, 0.);
	for ( size_t i = first; i <= last; ++i ){
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
		inverse_beta(args.ns->inversebeta) {
	magnetic_windC = vecd(Nx());
	const double C0 = inverse_beta * 2*mu_magn*mu_magn / GM;
	for (size_t i = 0; i < Nx_; i++) {
		double brackets;
		if (R_[i] < R_cor) {
			brackets = -7. + 4. * std::pow(R_[i] / R_cor, 1.5);
		} else {
			brackets = 7. - 10. * std::pow(R_cor / R_[i], 1.5);
		}
		magnetic_windC[i] += C0 / (R_[i]*R_[i]*R_[i]*R_[i]) * brackets;
	}
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
	for (ii = 0; ii < Nx() - 1; ii++) {
		if (R().at(ii+1) > R_m){
			break;
		}
	}
	if (ii >= Nx() - 2) {
		throw std::runtime_error("R_in > R_out");
	}
	R_m = R().at(ii);

	double new_F_in;
	if (R_m < R_cor) {
		const double n_ws = 1 - k_t * xi_pow_minus_7_2 * std::pow(R_m / R_cor, 3);
		new_F_in = (1 - n_ws) * Mdot_in() * std::sqrt(GM * R_m);
	} else {
		new_F_in = args_ns->Fdead * std::pow(R_dead / R_m, 3);
	}
	F_in_ = F_[0] = new_F_in;

	if (ii > 0) {
		Nx_ -= ii;
		h_.erase(h_.begin(), h_.begin() + ii);
		R_.erase(R_.begin(), R_.begin() + ii);
		F_.erase(F_.begin(), F_.begin() + ii);
	}
}

const vecd FreddiNeutronStarEvolution::windC() const {
	auto C = FreddiEvolution::windC();
	for (size_t i = 0; i < C.size(); i++) {
		C[i] += magnetic_windC[i];
	}
	return C;
}
