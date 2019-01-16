#include "freddi_state.hpp"

#include <cmath>
#include <exception>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "constants.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"


FreddiState::FreddiState(const FreddiArguments& args_, wunc_t wunc_):
		args(args_),
		Nt(static_cast<size_t>(std::round(args_.calc->time / args_.calc->tau))),
		GM(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args_.basic->Mx),
		eta(efficiency_of_accretion(args_.basic->kerr)),
		cosi(std::cos(args_.flux->inclination / 180 * M_PI)),
		cosiOverD2(std::cos(args_.flux->inclination / 180 * M_PI) / (args_.flux->distance * args_.flux->distance)),
		oprel(*args_.disk->oprel.get()),
		wunc(wunc_),
		Nx_(args_.calc->Nx),
		h_(Nx_),
		R_(Nx_),
		F_(Nx_),
		Mdot_out_(args.disk->Mdotout) {
	initializeGrid();
	initializeF();
}


void FreddiState::initializeGrid() {
	const double h_in = args.basic->h(args.basic->rin);
	const double h_out = args.basic->h(args.basic->rout);

	for (size_t i = 0; i < Nx_; ++i) {
		if (args.calc->gridscale == "log") {
			h_[i] = h_in * pow(h_out / h_in, i / (Nx_ - 1.));
		} else if (args.calc->gridscale == "linear") {
			h_[i] = h_in + (h_out - h_in) * i / (Nx_ - 1.);
		} else {
			throw std::logic_error("Wrong gridscale");
		}
		R_[i] = h_[i] * h_[i] / GM;
	}
}


void FreddiState::initializeF() {
	const double h_in = args.basic->h(args.basic->rin);
	const double h_out = args.basic->h(args.basic->rout);

	if (args.disk->initialcond == "power" or args.disk->initialcond == "powerF") {
		for (size_t i = 0; i < Nx_; ++i) {
			F_[i] = args.disk->F0 * std::pow((h_[i] - h_in) / (h_out - h_in), args.disk->powerorder);
		}
	} else if (args.disk->initialcond == "powerSigma") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double Sigma_to_Sigmaout = std::pow((h_[i] - h_in) / (h_out - h_in), args.disk->powerorder);
			F_[i] = args.disk->F0 * std::pow(h_[i] / h_out, (3. - oprel.n) / (1. - oprel.m)) *
					std::pow(Sigma_to_Sigmaout, 1. / (1. - oprel.m));
		}
	} else if (args.disk->initialcond == "sinusF" || args.disk->initialcond == "sinus") {
		for (size_t i = 0; i < Nx_; ++i) {
			F_[i] = args.disk->F0 * std::sin((h_[i] - h_in) / (h_out - h_in) * M_PI_2);
		}
	} else if (args.disk->initialcond == "quasistat") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double xi_LS2000 = h_[i] / h_out;
			F_[i] = args.disk->F0 * oprel.f_F(xi_LS2000) * (1. - h_in / h_[i]) / (1. - h_in / h_out);
		}
	} else if (args.disk->initialcond == "gaussF") {
		for (size_t i = 0; i < Nx_; ++i) {
			const double xi = (h_[i] - h_in) / (h_out - h_in);
			F_[i] = args.disk->F0 * exp(-(xi - args.disk->gaussmu) * (xi - args.disk->gaussmu) /
										 (2. * args.disk->gausssigma * args.disk->gausssigma));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}

	if (args.disk->wind == "no") {
		wind_.reset(new NoWind(this));
	} else if (args.disk->wind == "SS73C") {
		wind_.reset(new SS73CWind(this));
	} else if (args.disk->wind == "Cambier2013") { // Cambier & Smith 1303.6218
		wind_.reset(new Cambier2013Wind(this));
	} else if (args.disk->wind == "__testA__") {
		wind_.reset(new __testA__Wind(this));
	} else if (args.disk->wind == "__testB__") {
		wind_.reset(new __testB__Wind(this));
	} else if (args.disk->wind == "__testC__") {
		wind_.reset(new __testC__Wind(this));
	} else if (args.disk->wind == "__testC_q0_Shields1986__") {
		wind_.reset(new __testC_q0_Shields1986__(this));
	} else {
		throw std::logic_error("Wrong wind");
	}
}


void FreddiState::step(double tau) {
	set_Mdot_in_prev();
	invalidate_disk_optional_structure();
	i_t_ ++;
	t_ += tau;
	wind_->update();
}


double FreddiState::Mdot_in() const {
	return (F()[1] - F()[0]) / (h()[1] - h()[0]);
}


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
	if (!disk_str_.Lx_) {
		disk_str_.Lx_ = (Luminosity(R_, Tph_X(), args.flux->emin, args.flux->emax, 100)
				/ pow(args.flux->colourfactor, 4.));
	}
	return *disk_str_.Lx_;
}


const vecd& FreddiState::W() {
	if (!disk_str_.W_) {
		disk_str_.W_ = wunc(h_, F_, 0, Nx_ - 1);
	}
	return *disk_str_.W_;
}


const vecd& FreddiState::Sigma() {
	if (!disk_str_.Sigma_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = W()[i] * GM * GM / (4. * M_PI * pow(h()[i], 3.));
		}
		disk_str_.Sigma_ = std::move(x);
	}
	return *disk_str_.Sigma_;
}


const vecd& FreddiState::Tph() {
	if (!disk_str_.Tph_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = std::pow(std::pow(Tph_vis()[i], 4) + Qx()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		disk_str_.Tph_ = std::move(x);
	}
	return *disk_str_.Tph_;
}


const vecd& FreddiState::Tirr() {
	if (!disk_str_.Tirr_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = std::pow(Qx()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		disk_str_.Tirr_ = std::move(x);
	}
	return *disk_str_.Tirr_;
}


const vecd& FreddiState::Qx() {
	if (!disk_str_.Qx_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = (Cirr()[i] * eta * Mdot_in() * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT
					/ (4. * M_PI * R()[i] * R()[i]));
		}
		disk_str_.Qx_ = std::move(x);
	}
	return *disk_str_.Qx_;
}


const vecd& FreddiState::Cirr() {
	if (!disk_str_.Cirr_) {
		if (args.irr->irrfactortype == "const") {
			disk_str_.Cirr_ = vecd(Nx_, args.irr->Cirr);
		} else if (args.irr->irrfactortype == "square") {
			vecd x(Nx_);
			for (size_t i = 0; i < Nx_; i++) {
				x[i] = args.irr->Cirr * (Height()[i] / R()[i]) * (Height()[i] / R()[i]);
			}
			disk_str_.Cirr_ = std::move(x);
		} else {
			throw std::logic_error("Wrong irrfactor");
		}
	}
	return *disk_str_.Cirr_;
}


const vecd& FreddiState::Height() {
	if (!disk_str_.Height_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = oprel.Height(R()[i], F()[i]);
		}
		disk_str_.Height_ = std::move(x);
	}
	return *disk_str_.Height_;
}


const vecd& FreddiState::Tph_vis() {
	if (!disk_str_.Tph_vis_) {
		vecd x(Nx_);
		for (size_t i = 0; i < Nx_; i++) {
			x[i] = (GM * std::pow(h()[i], -1.75)
					* std::pow(3. / (8. * M_PI) * F()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25));
		}
		disk_str_.Tph_vis_ = std::move(x);
	}
	return *disk_str_.Tph_vis_;
}


const vecd& FreddiState::Tph_X() {
	if (!disk_str_.Tph_X_) {
		vecd x(Nx_);
		x[0] = 0;
		for (size_t i = 1; i < Nx_; i++) {
			x[i] = (args.flux->colourfactor
					* T_GR(R()[i], args.basic->kerr, args.basic->Mx, F()[i]
					/ (h()[i] - h()[0]), R()[0]));
		}
		disk_str_.Tph_X_ = std::move(x);
	}
	return *disk_str_.Tph_X_;
}


double FreddiState::lazy_magnitude(boost::optional<double>& m, double lambda, double F0) {
	if (!m){
		m = magnitude(lambda, F0);
	}
	return *m;
}


FreddiState::BasicWind::BasicWind(const FreddiState *state):
		state(state),
		A_(state->Nx(), 0.), B_(state->Nx(), 0.), C_(state->Nx(), 0.) {}

FreddiState::SS73CWind::SS73CWind(const FreddiState *state):
		BasicWind(state) {
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state->GM;
	const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT *
			state->eta);
	for (size_t i = 0; i < state->Nx_; ++i) {
		C_[i] = -Mdot_crit / (2 * M_PI * state->R_.front() * state->R_[i]) *
				(4 * M_PI * state->h_[i] * state->h_[i] * state->h_[i]) / (state->GM * state->GM);
	}
}

FreddiState::Cambier2013Wind::Cambier2013Wind(const FreddiState *state):
		BasicWind(state),
		kC(state->args.disk->windparams.at(0)),
		R_IC2out(state->args.disk->windparams.at(1)) {
	const auto disk = state->args.disk;
	const double m_ch0 = -kC * disk->Mdot0 / (M_PI * state->R_.back()*state->R_.back());  // dM / dA
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state->GM;
	const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * state->eta);
	const double eta = 0.025 * 33 * disk->Mdot0 / Mdot_crit;
	const double R_iC = R_IC2out * state->R_.back();
	for (size_t i = 0; i < state->Nx_; ++i) {
		const double xi = state->R_[i] / R_iC;
		const double C0 = m_ch0 * (4 * M_PI * state->h_[i]*state->h_[i]*state->h_[i]) / (state->GM*state->GM);
		C_[i] = C0 * std::pow((1 + std::pow(0.125 * eta / xi, 2))
								  / (1 + 1. / (std::pow(eta, 8) * std::pow(1 + 262 * xi * xi, 2))), 1. / 6.)
					* std::exp(-std::pow(1 - 1. / std::sqrt(1 + 0.25 / (xi * xi)), 2) / (2 * xi));
	}
}

FreddiState::__testA__Wind::__testA__Wind(const FreddiState *state):
		BasicWind(state),
		kA(state->args.disk->windparams.at(0)) {
	const double A0 = -kA / ((state->h_.back() - state->h_.front()) * (state->h_.back() - state->h_.front()));
	for (size_t i = 0; i < state->Nx_; ++i) {
		A_[i] = A0 * (state->h_[i] - state->h_.front());
	}
}

FreddiState::__testB__Wind::__testB__Wind(const FreddiState *state):
		BasicWind(state),
		kB(state->args.disk->windparams.at(0)) {
	const double B0 = -kB / ((state->h_.back() - state->h_.front()) * (state->h_.back() - state->h_.front()));
	for (size_t i = 0; i < state->Nx_; ++i) {
		B_[i] = B0;
	}
}

FreddiState::__testC__Wind::__testC__Wind(const FreddiState *state):
		BasicWind(state),
		kC(state->args.disk->windparams.at(0)) {
	const double C0 = kC * state->args.disk->Mdotout / (state->h_.back() - state->h_.front());
	const double h_wind_min = state->h_.back() / 2;
	for (size_t i = 0; i < state->Nx_; ++i) {
		if (state->h_[i] > h_wind_min) {
			C_[i] = C0 * 0.5 *
					(std::cos(2. * M_PI * (state->h_[i] - h_wind_min) / (state->h_.back() - h_wind_min)) - 1);
		}
	}
}

FreddiState::__testC_q0_Shields1986__::__testC_q0_Shields1986__(const FreddiState *state):
		BasicWind(state),
		kC(state->args.disk->windparams.at(0)),
		R_windmin2out(state->args.disk->windparams.at(1)) {}

void FreddiState::__testC_q0_Shields1986__::update() {
	BasicWind::update();
	const double h_wind_min = std::sqrt(R_windmin2out) * state->h_.back();
	for (size_t i = 0; i < state->Nx_; ++i) {
		if (state->h_[i] > h_wind_min) {
			C_[i] = -0.5/M_PI * kC * state->Mdot_in() /
					(std::log(1 / R_windmin2out) * state->R_[i] * state->R_[i]) *
					(4 * M_PI * state->h_[i]*state->h_[i]*state->h_[i]) / (state->GM*state->GM);
		}
	}
}
