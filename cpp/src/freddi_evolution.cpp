#include "freddi_evolution.hpp"

#include <cmath>
#include <string>

#include "arguments.hpp"
#include "nonlinear_diffusion.hpp"

using namespace std::placeholders;


FreddiEvolution::FreddiEvolution(const FreddiArguments &args):
		FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)) {}


void FreddiEvolution::step(const double tau) {
	truncateInnerRadius();
	FreddiState::step(tau);
	nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), Mdot_out(),
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last());
	truncateOuterRadius();
	star_.set_sources(star_irr_sources());
}


void FreddiEvolution::truncateOuterRadius() {
	if (args().disk->Thot <= 0. ){
		return;
	}
	// check proper value of accretion rate:
	if ((!std::isfinite(Mdot_in_prev())) || ( Mdot_in() < 0.0 ) || ( Mdot_in_prev() < 0.0 )) {
		return;
	}
	// check that Mdot decaying:
	if (Mdot_in() > Mdot_in_prev()) {
		return;
	}

	auto ii = last() + 1;
	if (Tirr().at(last()) / Tph_vis().at(last()) < args().disk->Tirr2Tvishot) {
	// when irradiation is not important
	// hot disc extends as far as Sigma>Sigma_max_cold(alpha_cold) and not farther than R_cooling_front and Tirr <= Thot 
		do {
			ii--;
			if (ii <= first()) throw std::runtime_error("Rout <= Rin");
		  } while( ( R().at(ii) > R_cooling_front ( R().at(ii)) ) && ( Sigma().at(ii) < Sigma_minus(R().at(ii)) ) && ( Tirr().at(ii) < args().disk->Thot ) );
		//} while( ( R().at(ii) > R_cooling_front ( R().at(ii)) ) && ( Sigma().at(ii) < Sigma_minus(R().at(ii)) ) );
	} else if (args().disk->boundcond == "Teff") {
	// irradiation is important, the boundary is at fixed Teff
		do {
			ii--;
			if (ii <= first()) throw std::runtime_error("Rout <= Rin");
		} while( Tph().at(ii) < args().disk->Thot );
	} else if (args().disk->boundcond == "Tirr") {
	// irradiation is important, the boundary is at fixed Tir
		do {
			ii--;
			if (ii <= first()) throw std::runtime_error("Rout <= Rin");
		} while( Tirr().at(ii) < args().disk->Thot );
	} else{
		throw std::invalid_argument("Wrong boundcond");
	}

	if ( ii <= last() - 1 ){
		current_.last = ii;
		invalidate_optional_structure();
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t _first, size_t _last) const {
	vecd W(_last + 1, 0.);
	for ( size_t i = _first; i <= _last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - oprel().m) * pow(h[i], oprel().n) / (1. - oprel().m) / oprel().D;
	}
	return W;
};
