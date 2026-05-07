#include "freddi_evolution.hpp"

#include <cmath>
#include <string>

#include "arguments.hpp"
#include "exceptions.hpp"
#include "nonlinear_diffusion.hpp"
#include <boost/math/tools/minima.hpp>

#define VERB_LEVEL_MESSAGES 30 

using namespace std::placeholders;


FreddiEvolution::FreddiEvolution(const FreddiArguments &args):
		FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)) {}


void FreddiEvolution::nonlinear_diffusion(const double tau) {
	nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), Mdot_outer_boundary(),
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last());
        
}

void FreddiEvolution::nonlinear_diffusion_outer_condition_depends_Mdotin(const double tau) {
    
    vecd Fvis_start = current_.F;

    //set_Mdot_outer_boundary(obtain_Mdot_outer_boundary());
    double current_Mdot_outer_boundary = Mdot_outer_boundary();
    //std::printf("(1) current_Mdot_outer_boundary: %e | Mdot_in: %e  %e\n",  current_Mdot_outer_boundary, Mdot_in(),current_Mdot_outer_boundary/Mdot_in() );
    
	std::uintmax_t maxit = 100;
	double left = 0.2;
	double right = 2.;
    
	std::cout.precision(10);
    int bits = 12;
    int iteration = 0;
    //double Mdot_out_from_F ;
    //double Mdot_in;
    // 2. Define the lambda (the function to minimize)
	auto objective_function = [this, &iteration, &Fvis_start, &current_Mdot_outer_boundary, tau](double coef) {
		iteration++;
		
        // Reset state for this specific guess
		current_.F = Fvis_start;
        //current_Mdot_outer_boundary = Mdot_outer_boundary();
        //std::printf("2a Mdot_in: %e Mdot_out: %e coef=%.6f\n",   Mdot_in(), current_Mdot_outer_boundary/Mdot_in(),coef);
        
        //printf (" Fs before: %e %e %e %e\n", current_.F[first()], current_.F[first()+1], current_.F[last()-1], current_.F[last()]);
		// Run the physics
		nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), coef*current_Mdot_outer_boundary,
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last()
		);


        //printf (" Fs  after: %e %e %e %e\n", current_.F[first()], current_.F[first()+1], current_.F[last()-1], current_.F[last()]);

        //double coef_after = current_Mdot_outer_boundary/Mdot_in();
        double Mdot_calc_from_F = Mdot_out_from_F();
        //Mdot_out_from_F = Mdot_calc_from_F;
        //Mdot_in = Mdot_in();
        //std::printf("2b Mdot_in: %e Mdot_out_calc: %e coef_after %e\n",   Mdot_in(),  Mdot_calc_from_F, coef_after);
        
		// Calculate from Mdot_in  appropriate Mdot_out_boundary
        set_Mdot_outer_boundary(obtain_Mdot_outer_boundary());

        double iter_Mdot_outer_boundary = Mdot_outer_boundary();
        

		double Answer = iter_Mdot_outer_boundary;
		// compare two values of Mdot_out_boundary
		double Answer_div = 1.0 - iter_Mdot_outer_boundary / Mdot_calc_from_F;
		
		//double result = std::pow(Answer_div, 2);
        double result = std::fabs(Answer_div);
		// 3. Print INSIDE the lambda so it logs every step Brent takes
		//std::printf("Iter: %3d| coef: %f | now_Mdot_outer_boundary: %e | result: %e | previous_Mdot_outer_bound: %e Mdot_out_from_F %e check %e\n", iteration, coef, Answer, result, current_Mdot_outer_boundary, Mdot_calc_from_F, Mdot_calc_from_F/Mdot_in());
        //std::printf("Iter: %3d| coef: %f \n", iteration, coef);
		return result; 
	};
    // 4. Now call the solver using the name of the lambda
	std::pair<double, double> brent_result = boost::math::tools::brent_find_minima(objective_function, left, right, bits, maxit);
    // The final result found by Brent
    
	//std::cout << "---------------SEARCH FINISHED. Best : " << brent_result.first << " Min Error: " << brent_result.second << std::endl;
    
    double m1 = Mdot_out_from_F();
    double m2 = Mdot_in();
    
    if (std::fabs((std::fabs(m1/m2) - args().disk->DIM_front_Mdot_factor)> 0.01 ))  {

         std::printf("(W) RESULT factor current_Mdot_outer_boundary/Mdot_in:  %e for Mdot_in=%e; required=%.3e\n", m1/m2,  m2, args().disk->DIM_front_Mdot_factor);
         //std::getchar();  
    }
    
}
		
void FreddiEvolution::step(const double tau) {
         
	//if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "cA_ t="<< [freddi]() {return sToDay(current_.t);}  <<"\n" << std::endl;}
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A_ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
	truncateInnerRadius();
	
	FreddiState::step(tau);
	
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A__ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
	//nonlinear_diffusion(tau);
	if (Mdot_outer_boundary() == 0.0) {
        nonlinear_diffusion(tau); // valid for outer boundary condition Mdotin=0
    } else {
        nonlinear_diffusion_outer_condition_depends_Mdotin(tau); //  iterations until the outer boundary condition fullfulls at the end of the step
    }
	
    /*nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), Mdot_outer_boundary(),
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last());*/
	
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_A___ t="<< sToDay(current_.t)  <<"\n" << std::endl;}
	truncateOuterRadius();
	//truncateInnerRadius();//added @@
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

bool FreddiEvolution::check_ring_is_cold(const int ii) {
    // returns true for cold and false for hot

    // R().at(ii), the radius where Mdot = 0, or ?
    // multiplied by radius_popravka, gives the hot zone radius
    // Sigma_minus is the maximum density on the cold branch
    // Sigma_plus is the minimum density on the hot branch     = Sigma_min(Menou+1999)


    if (args().disk->boundcond == "Teff") {
        if (Tph().at(ii) >= args().disk->Thot) {
            return false; // HOT
        } else {
            return true; // COLD
        }
    }

    // determine radius, inside which irradiation is not significant:
    find_R_max_where_Qirr_works();
    
    if ((args().disk->boundcond == "Tirr") || (args().disk->boundcond == "Tirr_Ham??eury")) {
        if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "boundcond=Tirr="<<Tirr_critical (R().at(ii), ii) << std::endl;}
        if (Tirr().at(ii) >= Tirr_critical (R().at(ii), ii)) {
            set_R_dotM0_before_shift( R().at(ii) );
            return false; // HOT
        } else {
            return true; // COLD

        }
    }

    if (!(args().disk->boundcond == "DIM")) {
        throw std::invalid_argument("Wrong boundcond at check_ring_is_cold() in freddi_state");
    }


    if (  Tirr().at(ii) >= Tirr_critical (R().at(ii), ii) ) {
        // if irradiation temperature is greater than critical, disc cannot be cold
        // If there is no scattering, the disc beyond r, where dotM=0, is in the shadow from the direct central radiation
        
        if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_K Tirr high - HOT \n" << std::endl;}
        // memorise radius before shift of the hot boundary
        set_R_dotM0_before_shift( R().at(ii) );
        return false; // HOT
    } else  {
        // check_Sigma_approach -> "cold-front-approach" ?
        // Menou99a -> Sigma_Menou99 ?
        // if (args().disk->check_Sigma_approach == "Menou99a") { ?
        // when irradiation stops preventing hot-zone shrinking,
        // then conditions at Rfront matter:
        
        // newest approach: radius_popravka = 1
        double radius_popravka = Rfront_Rhot( R().at(ii), oprel().Height(R()[ii], F()[ii])/R().at(ii));


        // calculate sigma_popravka:
        double sigma_factor;
        if (radius_popravka == 1. ) {
            // check the conditions at r where Mdot=0
            sigma_factor = 1.;
        } else {
            if (args().disk->check_Sigma_approach == "Menou99a") {
                sigma_factor = 1./4.3 ; // See fig.8 of Menou+1999; this correctly work only for constant radius_popravka
            } else if (args().disk->check_Sigma_approach == "simple") {
                sigma_factor = pow(radius_popravka, -0.75);
            } else if (args().disk->check_Sigma_approach == "critical_Teff") {

            } else {
                throw std::invalid_argument("Wrong check_Sigma_approach [Menou99a/simple/critical_Teff]");
            }
        }

        // check if surface density at front is larger than critical cold Sigma_minus
        //
        if (Sigma().at(ii)/pow(radius_popravka, -0.75) >  Sigma_minus(R().at(ii) * radius_popravka)) {
            // disc cannot be cold if density is greater than critical for cold state
            if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_L pop ="<<  radius_popravka <<  " Sigma is high - HOT \n" << std::endl;}
            set_R_dotM0_before_shift( R().at(ii) );
            return false; // HOT
        } else {

            if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_S Sigma is  LOW -> cooling wave OR collapse ii="<<ii<<" pop =" << radius_popravka <<"sigma_minus=" << Sigma_minus(R().at(ii) * radius_popravka) <<  std::endl;}

            // check if surface density at front is less than critical hot Sigma_plus:
            if ( Sigma().at(ii)*sigma_factor < Sigma_plus(R().at(ii)*radius_popravka) ) {
                //ring cannot be hot if density is lower than critical for hot state
                // and if cooling front had time to reach this radius
                if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_M Sigma is < Sigma_hot_crit sigma_outer ="<<  last() << " -- " << args().calc->Nx-1 <<"\n" << std::endl;}
                // last = Nx-1 means that code collapses to a starting Rhot : all outer disc is cast COLD
                if ( (last() == args().calc->Nx-1) || (current_.R_dotM0_before_shift == 0.) ) {
                    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_O pop ="<<  radius_popravka << " Sigma is low and R_dotM0_before_shift=" << current_.R_dotM0_before_shift<< "- COLD \n" << std::endl;}
//                     return 1; // COLD
                } else {
                    //throw std::invalid_argument("check_ring_is_cold: logic mistake 1");
                }
                return 1; // COLD
            }
// 		if ( last() < args().calc->Nx-1) {
// 		    if (Sigma_plus(R().at(ii+1))== 0.)  { Sigma().at(ii)
// 			// Sigma is low BUT ring is HOT
// 			// that means that the front is propagating and it is propagating with
// 			// meaningful velocity.
// 			// otherwise we search for a starting Rhot which position can be anything
// 			if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_N pop ="<<  radius_popravka << " Sigma is low BUT HOT  " <<  Sigma_plus(R().at(ii+1)) << "\n"<< std::endl;}
// 			return 0;
// 		    }
// 		}
// 		r????eturn 1;

            //  R_cooling_front = r - v_cooling_front(r, sigma_at_r) * args().calc->tau;
//   		if (radius_popravka <= args().disk->Rfront_Mdotzero_factor) {
// 		    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_PP pop ="<< radius_popravka <<  " radius_popravka < Rfront_Mdotzero_factor->COLD \n" << std::endl;}
//   		    return 1;
//   		}

            
            // here below Sigma > Sigma_plus
            
            
            if  (current_.R_dotM0_before_shift == 0 ) {   // and 
                // COLD ; just collapse to initial position:
                if (last() != args().calc->Nx-1) {
                    throw std::invalid_argument("check_ring_is_cold: logic mistake 1.5");
                }
                return 1;  // COLD
            }
            
            
            if ((radius_popravka * R().at(ii) > R_cooling_front ( Rfront_Rhot(R()[last()], oprel().Height(R()[last()], F()[last()])/R()[last()] ) * R()[last()] , Sigma().at(last())*sigma_factor) )  && ( last() < args().calc->Nx-1)  ){
                //radius is beyond front
                if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "++c_P pop ="<< radius_popravka <<  " beyond front - COLD \n" << std::endl;}
                return 1;// COLD
            } else {
                if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "++c_Q pop ="<<  radius_popravka <<  " inwards front - HOT \n" << std::endl;}
                // in fact front starts from a greater distance?
                set_R_dotM0_before_shift( R().at(ii) );
                return false;// HOT
            }
        }
    }
    throw std::invalid_argument("check_ring_is_cold: logic mistake 2");
    return false;
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t _first, size_t _last) const {
	vecd W(_last + 1, 0.);
	for ( size_t i = _first; i <= _last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - oprel().m) * pow(h[i], oprel().n) / (1. - oprel().m) / oprel().D;
	}
	return W;
};


#undef VERB_LEVEL_MESSAGES
