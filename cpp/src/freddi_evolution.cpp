#include "freddi_evolution.hpp"

#include <cmath>
#include <string>

#include "arguments.hpp"
#include "exceptions.hpp"
#include "nonlinear_diffusion.hpp"

#define VERB_LEVEL_MESSAGES 30 

using namespace std::placeholders;


FreddiEvolution::FreddiEvolution(const FreddiArguments &args):
		FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)) {}


void FreddiEvolution::step(const double tau) {
	//if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "cA_ t="<< [freddi]() {return sToDay(current_.t);}  <<"\n" << std::endl;}
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A_ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
	truncateInnerRadius();
	FreddiState::step(tau);
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A__ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
	nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), Mdot_outer_boundary(),
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last());
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A___ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
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
	
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout  <<"c_A last ii =  "<< last()  <<" last R=" <<R()[last()] <<" \n" << std::endl;}
	
	
	do {
	    ii--;
	    if (ii <= first()) throw RadiusCollapseException();
	}  while ( check_ring_is_cold(ii) );
	

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


#undef VERB_LEVEL_MESSAGES
