#include "freddi_evolution.hpp"

#include <cmath>
#include <string>

#include "arguments.hpp"
#include "exceptions.hpp"
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

int FreddiState::ring_state_vertical(const int ii) {
    // returns 1 for hot, 0 for cold
    
    double radius_popravka = args().disk->Rhot_Mdotzero_factor;
    
    // R().at(ii) is the radius where Mdot = 0
    // multiplied by radius_popravka, gives the hot zone radius
    // Sigma_minus is the maximum density on the cold branch  
    // Sigma_plus is the minimum density on the hot branch     = Sigma_min(Menou+1999)
    
    if ( Tirr().at(ii)/pow(radius_popravka,0.5) > args().disk->Thot) {
	// if irradiation temperature is greater than critical, disc cannot be cold
	return 1;
    } else {
	// check if surface density at front is larger than critical cold Sigma_minus
	// 
	if (Sigma().at(ii)/pow(radius_popravka, -0.75) >  Sigma_minus(R().at(ii) * radius_popravka)) {
	    // disc cannot be cold if density is greater than critical for cold state
	    return 1;
	} else {
	    // check if surface density at front is less than critical hot Sigma_plus
	    double sigma_factor;
	    if (args().disk->check_Sigma_approach == "Menou99a") {
		sigma_factor = 4.3 ; // See fig.8 of Menou+1999 
	    } else if (args().disk->check_Sigma_approach == "simple") {
		sigma_factor = pow(radius_popravka, -0.75);
	    } else {
		throw std::invalid_argument("Wrong check_Sigma_approach [Menou99a/simple]");
	    }
	    if ( Sigma().at(ii)/sigma_factor < Sigma_plus(R().at(ii)*radius_popravka) ) {
		//disc cannot be hot if density is lower than critical for hot state
		return 0;
	    } else {
		// check cooling front position
		if (radius_popravka * R().at(ii) > R_cooling_front ( radius_popravka*R().at(ii)) ) {
		    //radius is beyond front 
		    return 0;
		} else {
		    return 1;
		}
	    }
	}
    }
    throw std::invalid_argument("ring_state_vertical: logic mistake");
    return 1;
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
	
	
	double radius_popravka =  args().disk->Rhot_Mdotzero_factor; //1.8; 2.1; 1.7;
	
	if (Tirr().at(last()) / Tph_vis().at(last()) < args().disk->Tirr2Tvishot/pow(radius_popravka,0.25)) {
	// when irradiation is not important
	// (A) hot disc extends as far as Sigma>Sigma_max_cold(alpha_cold) and not farther than R_cooling_front and Tirr >= Thot  and Teff_vis >=Teff_plus (condition below is the opposite)
	// or (B) hot disk cannot exist for Sigma< Sigma_crit_hot (=Sigma+) while Tirr > Tcrit: this condition ensures fast
	// convergence to small disc radius in the case without irradiation    
	    //DEBUG std::cout << Tirr().at(last()) << "\n" << std::endl;
	    
	    
		do {
			ii--;
			if (ii <= first()) throw RadiusCollapseException();
		  } while( 
		        // conditions for cold zone:
		        //(A) (1)  hot disc extends not farther than R_cooling_front
		        // sent to DIM team
		        ((args().disk->check_state_approach == "before2024") &&
				(
				    ( radius_popravka * R().at(ii) > R_cooling_front ( radius_popravka*R().at(ii)) ) 
				    
				// Sigma_minus is the maximum density on the cold branch  
				// Sigma_plus is the minimum density on the hot branch     = Sigma_min(Menou+1999)
				// (2) hot disc extends as far as Sigma>Sigma_max_cold(alpha_cold)    
	// 			&& ( Sigma().at(ii)/pow(radius_popravka,0.75) < Sigma_minus(R().at(ii)) ) 
			//	&& ( Sigma().at(ii) < Sigma_minus(R().at(ii)) ) 
			//	&& ( Sigma().at(ii)/4.3 < Sigma_minus(R().at(ii)) *  pow(radius_popravka,1.18) )  // error?
				&& ((args().disk->check_Sigma_approach == "Menou99a") && ( Sigma().at(ii)/4.3 < Sigma_plus(R().at(ii)) *  pow(radius_popravka,1.11) ))  // see Fig. 8 of Menou+99
				
	//			 && ( Tirr().at(ii)/pow(radius_popravka,0.5) < args().disk->Thot ) 
				// && ( Tirr().at(ii)/pow(radius_popravka,0.5) < 9040. - 2216.* 1./(pow(Tirr().at(ii) / Tph_vis().at(ii),4.)*radius_popravka) )   // this is valid only if Tirr>Tvis
				
				// (3) in the cold disc temperature < critical
				&& ( Tirr().at(ii)/pow(radius_popravka,0.5) < args().disk->Thot) //6800 was  
				)) 
			|| 
			( (args().disk->check_state_approach == "logic") && ( ring_state_vertical(ii) == 0) )
			
			// (B)
			// || (
			//    ( Sigma().at(ii) < Sigma_plus(R().at(ii)) ) 
			//    && ( Tirr().at(ii) < args().disk->Thot ) 
			//)
			
			// (3) next line is added according to Tavleev+23 results: it overrides the line above if Thot=10000. In fact, the line above should be deleted.
		//	&& ( Tirr().at(ii) < 9040. - 2216. ) 
			// (4 - possibly wrong) next line is added to prevent cooling wave when effective viscous temperature is too high:
			//&& ( Tph_vis().at(ii) < Teff_plus(R().at(ii)) )
		);
		//} while( ( R().at(ii) > R_cooling_front ( R().at(ii)) ) && ( Sigma().at(ii) < Sigma_minus(R().at(ii)) ) );
	} else if (args().disk->boundcond == "Teff") {
	// irradiation is important, the boundary is at fixed Teff
		do {
			ii--;
			if (ii <= first()) throw RadiusCollapseException();
		} while( Tph().at(ii) < args().disk->Thot );
		    
		    
	} else if (args().disk->boundcond == "Tirr") {
	    
	    double Rcheck;
	    Rcheck = std::min( R().at(ii-1)*radius_popravka, args().basic->rout);
	    double radius_popravka_variable = Rcheck/R().at(ii-1);
	    // DEBUG std::cout << radius_popravka_variable << " " << 9040. - 2216.* 1./(pow(Tirr().at(ii-1) / Tph_vis().at(ii-1),4.)*radius_popravka_variable) <<  "\n" << std::endl;
	// irradiation is important, the boundary is at fixed Tir
		do {
			ii--;
			if (ii <= first()) throw RadiusCollapseException();
		//} while( Tirr().at(ii) < 9040. - 2216.* Tph_vis().at(ii) / Tirr().at(ii) ); 
		    // according to Tavleev+23 results; it is applicable in quasi-stationary disc (without cooling wave)
// 		    std::cout << Tirr().at(ii) << "\n" << std::endl;
// 		    std::cout << args().disk->Thot << "\n" << std::endl;
		// old: 
// 		 } while( Tirr().at(ii)/pow(radius_popravka,0.5) < args().disk->Thot);
		 } while( Tirr().at(ii)/pow(radius_popravka_variable,0.5) < args().disk->Thot);
//		} while( Tirr().at(ii)/pow(radius_popravka_variable,0.5) < 9040. - 2216.* 1./(pow(Tirr().at(ii) / Tph_vis().at(ii),4.)*radius_popravka_variable) );  //variable Tcrit does not correspond to SIM model
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
