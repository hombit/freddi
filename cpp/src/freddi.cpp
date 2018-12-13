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


FreddiEvolution::FreddiEvolution(const FreddiArguments &args_):
		Nt(static_cast<size_t>(std::round(args_.calc->time / args_.calc->tau))),
		args(&args_),
		GM(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args_.basic->Mx),
		eta(efficiency_of_accretion(args_.basic->kerr)),
		cosi(std::cos(args_.flux->inclination / 180 * M_PI)),
		cosiOverD2(std::cos(args_.flux->inclination / 180 * M_PI) / (args_.flux->distance * args_.flux->distance)),
		oprel(args_.disk->oprel.get()),
		wunc(std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)),
		state_(new FreddiState(this)) {}


void FreddiEvolution::set_Mdot_in_prev() {
	set_Mdot_in_prev(state_->Mdot_in());
}


void FreddiEvolution::step(const double tau) {
	set_Mdot_in_prev();
	updateWind();
	state_.reset(new FreddiState(*state_, tau));
	nonlinear_diffusion_nonuniform_wind_1_2(
			args->calc->tau, args->calc->eps,
			state_->F_in(), state_->Mdot_out(),
			state_->windA(), state_->windB(), state_->windC(),
			wunc,
			state_->h(), state_->F_);
	truncateOuterRadius();
}


std::vector<FreddiState> FreddiEvolution::evolve() {
	std::vector<FreddiState> states;
	states.push_back(*state_);
	for (size_t i_t = 0; i_t < Nt; i_t++) {
		step();
		states.push_back(*state_);
	}
	return states;
}


void FreddiEvolution::updateWind() {
	if (args->disk->wind == "__testC_q0_Shields1986__") {
		const double kC = 1;
		const double r_dimless_wind_min = 0.9;
		const double h_wind_min = std::sqrt(r_dimless_wind_min) * state_->h_.back();
		for (size_t i = 0; i < state_->Nx_; ++i) {
			if (state_->h_[i] > h_wind_min) {
				state_->windC_[i] = -0.5/M_PI * kC * state_->Mdot_in()
						/ (std::log(1 / r_dimless_wind_min) * state_->R_[i] * state_->R_[i])
						* (4 * M_PI * state_->h_[i]*state_->h_[i]*state_->h_[i]) / (GM*GM);
			}
		}
	}
}


void FreddiEvolution::truncateOuterRadius() {
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
	if (state_->Mdot_in() > Mdot_in_prev) {
		return;
	}

	auto ii = state_->Nx_;
	if (args->disk->boundcond == "Teff") {
		do {
			ii--;
		} while( state_->Tph().at(ii) < args->disk->Thot );
	} else if (args->disk->boundcond == "Tirr") {
		do {
			ii--;
		} while( state_->Tirr().at(ii) < args->disk->Thot );
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii < state_->Nx_-1 ){
		state_->Nx_ = ii+1;
		// F.at(Nx-2) = F.at(Nx-1) - Mdot_out / (2.*M_PI) * (h.at(Nx-1) - h.at(Nx-2));
		state_->h_.resize(state_->Nx_);
		state_->R_.resize(state_->Nx_);
		state_->F_.resize(state_->Nx_);
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t first, size_t last) const {
	vecd W(last + 1, 0.);
	for ( size_t i = first; i <= last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - oprel->m) * pow(h[i], oprel->n) / (1. - oprel->m) / oprel->D;
	}
	return W;
};


// Equation from Lasota, Dubus, Kruk A&A 2008, Menou et al. 1999. Sigma_cr is from their fig 8 and connected to point where Mdot is minimal.
double FreddiEvolution::Sigma_hot_disk(double r) const {
	return 39.9 * pow(args->basic->alpha/0.1, -0.80) * pow(r/1e10, 1.11) * pow(args->basic->Mx/GSL_CONST_CGSM_SOLAR_MASS, -0.37);
};



FreddiNeutronStarEvolution::FreddiNeutronStarEvolution(const FreddiNeutronStarArguments &args):
		FreddiEvolution(args),
		args_ns(args.ns.get()),
		xi_pow_minus_7_2(std::pow(xi, -3.5)),
		R_m_min(std::max(args.ns->Rx, args.basic->rin)),
		mu_magn(0.5 * args.ns->Bx * args.ns->Rx*args.ns->Rx*args.ns->Rx),
		R_dead(std::cbrt(mu_magn*mu_magn / args.ns->Fdead)),
		R_cor(std::cbrt(GM / (4 * M_PI*M_PI * args.ns->freqx*args.ns->freqx))) {}


void FreddiNeutronStarEvolution::step(double tau) {
	set_Mdot_in_prev();
	state_.reset(new FreddiState(*state_, tau));
	truncateInnerRadius();
	nonlinear_diffusion_nonuniform_wind_1_2(
			args->calc->tau, args->calc->eps,
			state_->F_in(), state_->Mdot_out(),
			state_->windA(), state_->windB(), state_->windC(),
			wunc,
			state_->h(), state_->F_
	);
	truncateOuterRadius();
}


void FreddiNeutronStarEvolution::truncateInnerRadius() {
	if (args_ns->Fdead <= 0.) {
		return;
	}
	if ( state_->Mdot_in() > Mdot_in_prev ) {
		return;
	}

	const double R_alfven = args_ns->epsilonAlfven *
			std::pow(mu_magn*mu_magn*mu_magn*mu_magn / (state_->Mdot_in()*state_->Mdot_in() * GM), 1./7.);
	double R_m = std::max(R_m_min, R_alfven);
	R_m = std::min(R_m, R_dead);
	size_t ii;
	for (ii = 0; ii < state_->Nx() - 1; ii++) {
		if (state_->R().at(ii+1) > R_m){
			break;
		}
	}
	if (ii >= state_->Nx() - 2) {
		throw std::runtime_error("R_in > R_out");
	}
	R_m = state_->R().at(ii);

	double new_F_in;
	if (R_m < R_cor) {
		const double n_ws = 1 - k_t * xi_pow_minus_7_2 * std::pow(R_m / R_cor, 3);
		new_F_in = (1 - n_ws) * state_->Mdot_in() * std::sqrt(GM * R_m);
	} else {
		new_F_in = args_ns->Fdead * std::pow(R_dead / R_m, 3);
	}
	state_->F_in_ = state_->F_[0] = new_F_in;

	if (ii > 0) {
		state_->Nx_ -= ii;
		state_->h_.erase(state_->h_.begin(), state_->h_.begin() + ii);
		state_->R_.erase(state_->R_.begin(), state_->R_.begin() + ii);
		state_->F_.erase(state_->F_.begin(), state_->F_.begin() + ii);
	}
}



FreddiState::FreddiState(const FreddiEvolution* freddi):
		freddi(freddi),
		Nx_(freddi->args->calc->Nx),
		h_(Nx_),
		R_(Nx_),
		F_(Nx_),
		windA_(Nx_, 0.),
		windB_(Nx_, 0.),
		windC_(Nx_, 0.),
		Mdot_out_(freddi->args->disk->Mdotout) {
	initializeGrid();
	initializeF();
	initializeWind();
}


void FreddiState::initializeGrid() {
	const auto args = freddi->args;
	const double h_in = args->basic->h(args->basic->rin);
	const double h_out = args->basic->h(args->basic->rout);

	for (size_t i = 0; i < Nx_; ++i) {
		if (args->calc->gridscale == "log") {
			h_[i] = h_in * pow(h_out / h_in, i / (Nx_ - 1.));
		} else if (args->calc->gridscale == "linear") {
			h_[i] = h_in + (h_out - h_in) * i / (Nx_ - 1.);
		} else {
			throw std::logic_error("Wrong gridscale");
		}
		R_[i] = h_[i] * h_[i] / freddi->GM;
	}
}


void FreddiState::initializeF() {
	const auto args = freddi->args;
	const auto oprel = freddi->oprel;

	const double h_in = args->basic->h(args->basic->rin);
	const double h_out = args->basic->h(args->basic->rout);

	if (args->disk->initialcond == "power" or args->disk->initialcond == "powerF") {
		for (size_t i = 0; i < Nx_; ++i) {
			F_[i] = args->disk->F0 * std::pow((h_[i] - h_in) / (h_out - h_in), args->disk->powerorder);
		}
	} else if (args->disk->initialcond == "powerSigma") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double Sigma_to_Sigmaout = std::pow((h_[i] - h_in) / (h_out - h_in), args->disk->powerorder);
			F_[i] = args->disk->F0 * std::pow(h_[i] / h_out, (3. - oprel->n) / (1. - oprel->m)) *
					std::pow(Sigma_to_Sigmaout, 1. / (1. - oprel->m));
		}
	} else if (args->disk->initialcond == "sinusF" || args->disk->initialcond == "sinus") {
		for (size_t i = 0; i < Nx_; ++i) {
			F_[i] = args->disk->F0 * std::sin((h_[i] - h_in) / (h_out - h_in) * M_PI_2);
		}
	} else if (args->disk->initialcond == "quasistat") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double xi_LS2000 = h_[i] / h_out;
			F_[i] = args->disk->F0 * oprel->f_F(xi_LS2000) * (1. - h_in / h_[i]) / (1. - h_in / h_out);
		}
	} else if (args->disk->initialcond == "gaussF") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double xi = (h_[i] - h_in) / (h_out - h_in);
			F_[i] = args->disk->F0 * exp(-(xi - args->disk->gaussmu) * (xi - args->disk->gaussmu) /
										 (2. * args->disk->gausssigma * args->disk->gausssigma));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}
}


void FreddiState::initializeWind() {
	const auto args = freddi->args;

	const double h_in = args->basic->h(args->basic->rin);
	const double h_out = args->basic->h(args->basic->rout);

	if (args->disk->wind == "no") {
		// Nothing to do here
	} else if (args->disk->wind == "SS73C") {
		const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
							 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * freddi->GM;
		const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * freddi->eta);
		for (size_t i = 0; i < Nx_; ++i) {
			windC_[i] = -Mdot_crit / (2 * M_PI * R_.front() * R_[i]) * (4 * M_PI * h_[i] * h_[i] * h_[i]) /
						(freddi->GM * freddi->GM);
		}
	} else if (args->disk->wind == "Cambier2013") { // Cambier & Smith 1303.6218
		const double kC = 3;
		const double m_ch0 = -kC * args->disk->Mdot0 / (M_PI * R_.back()*R_.back());  // dM / dA
		const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
							 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * freddi->GM;
		const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * freddi->eta);
		const double eta = 0.025 * 33 * args->disk->Mdot0 / Mdot_crit;
		const double R_iC = 1. * R_.back();
		for (size_t i = 0; i < Nx_; ++i) {
			const double xi = R_[i] / R_iC;
			const double C0 = m_ch0 * (4 * M_PI * h_[i]*h_[i]*h_[i]) / (freddi->GM*freddi->GM);
			windC_[i] = C0 * std::pow((1 + std::pow(0.125 * eta / xi, 2))
									  / (1 + 1. / (std::pow(eta, 8) * std::pow(1 + 262 * xi * xi, 2))), 1. / 6.)
						* std::exp(-std::pow(1 - 1. / std::sqrt(1 + 0.25 / (xi * xi)), 2) / (2 * xi));
		}
	} else if (args->disk->wind == "__testA__") {
		const double kA = 10;
		const double A0 = -kA / ((h_out - h_in) * (h_out - h_in));
		for (size_t i = 0; i < Nx_; ++i) {
			windA_[i] = A0 * (h_[i] - h_in);
		}
	} else if (args->disk->wind == "__testB__") {
		const double kB = 16.0;
		const double B0 = -kB / ((h_out - h_in) * (h_out - h_in));
		for (size_t i = 0; i < Nx_; ++i) {
			windB_[i] = B0;
		}
	} else if (args->disk->wind == "__testC__") {
		const double kC = 3.0;
		const double C0 = kC * args->disk->Mdotout / (h_out - h_in);
		const double h_wind_min = h_out / 2;
		for (size_t i = 0; i < Nx_; ++i) {
			if (h_[i] > h_wind_min) {
				windC_[i] = C0 * 0.5 * (std::cos(2. * M_PI * (h_[i] - h_wind_min) / (h_out - h_wind_min)) - 1);
			}
		}
	} else if (args->disk->wind == "__testC_q0_Shields1986__") {
		// To be set in updateWind()
	} else {
		throw std::logic_error("Wrong wind");
	}
}


FreddiState::FreddiState(const FreddiState& other, const double tau):
		freddi(other.freddi),
		Nx_(other.Nx_),
		h_(other.h_),
		R_(other.R_),
		F_(other.F_),
		t_(other.t_ + tau),
		i_t_(other.i_t_ + 1),
		Mdot_out_(other.Mdot_out_),
		windA_(other.windA_), windB_(other.windB_), windC_(other.windC_) {}


double FreddiState::integrate(const vecd& values) const {
	double integral = 0.;
	const size_t N = std::min(values.size(), R_.size());
	double stepR;
	for ( int i = 0; i < N; ++i ){
		if ( i == 0            ) stepR = R_[i+1] - R_[i  ];
		if ( i == N-1          ) stepR = R_[i  ] - R_[i-1];
		if ( i > 1 and i < N-1 ) stepR = R_[i+1] - R_[i-1];
		integral += 0.5 * values[i] * 2.*M_PI * R_[i] * stepR;
	}
	return integral;
}


double FreddiState::lazy_integrate(boost::optional<double>& x, const FreddiState::vecd& values) {
	if (!x) {
		x = integrate(values);
	}
	return *x;
}


double FreddiState::Lx() {
	if (!Lx_) {
		Lx_ = (Luminosity(R_, Tph_X(), freddi->args->flux->emin, freddi->args->flux->emax, 100)
				/ pow(freddi->args->flux->colourfactor, 4.));
	}
	return *Lx_;
}


const vecd& FreddiState::W() {
	if (!W_) {
		W_ = freddi->wunc(h_, F_, 0, Nx_ - 1);
	}
	return *W_;
}


const vecd& FreddiState::Sigma() {
	if (!Sigma_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = W()[i] * freddi->GM * freddi->GM / (4. * M_PI * pow(h()[i], 3.));
		}
		Sigma_ = std::move(x);
	}
	return *Sigma_;
}


const vecd& FreddiState::Tph() {
	if (!Tph_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = std::pow(std::pow(Tph_vis()[i], 4) + Qx()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		Tph_ = std::move(x);
	}
	return *Tph_;
}


const vecd& FreddiState::Tirr() {
	if (!Tirr_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = std::pow(Qx()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		Tirr_ = std::move(x);
	}
	return *Tirr_;
}


const vecd& FreddiState::Qx() {
	if (!Qx_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = (Cirr()[i] * freddi->eta * Mdot_in() * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT
					/ (4. * M_PI * R()[i] * R()[i]));
		}
		Qx_ = std::move(x);
	}
	return *Qx_;
}


const vecd& FreddiState::Cirr() {
	if (!Cirr_) {
		if (freddi->args->irr->irrfactortype == "const") {
			Cirr_ = vecd(Nx_, freddi->args->irr->Cirr);
		} else if (freddi->args->irr->irrfactortype == "square") {
			vecd x(Nx_);
			for (size_t i = 0; i < Nx_; i++) {
				x[i] = freddi->args->irr->Cirr * (Height()[i] / R()[i]) * (Height()[i] / R()[i]);
			}
			Cirr_ = std::move(x);
		} else {
			throw std::logic_error("Wrong irrfactor");
		}
	}
	return *Cirr_;
}


const vecd& FreddiState::Height() {
	if (!Height_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = freddi->oprel->Height(R()[i], F()[i]);
		}
		Height_ = std::move(x);
	}
	return *Height_;
}


const vecd& FreddiState::Tph_vis() {
	if (!Tph_vis_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = (freddi->GM * std::pow(h()[i], -1.75)
					* std::pow(3. / (8. * M_PI) * F()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25));
		}
		Tph_vis_ = std::move(x);
	}
	return *Tph_vis_;
}


const vecd& FreddiState::Tph_X() {
	if (!Tph_X_) {
		vecd x(Nx_);
		x[0] = 0;
		for (size_t i = 1; i < Nx_; i++) {
			x[i] = (freddi->args->flux->colourfactor
					* T_GR(R()[i], freddi->args->basic->kerr, freddi->args->basic->Mx, F()[i]
					/ (h()[i] - h()[0]), R()[0]));
		}
		Tph_X_ = std::move(x);
	}
	return *Tph_X_;
}


double FreddiState::lazy_magnitude(boost::optional<double>& m, double lambda, double F0) {
	if (!m){
		m = magnitude(lambda, F0);
	}
	return *m;
}
