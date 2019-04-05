#include "freddi_evolution.hpp"

#include <cmath>
#include <exception>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "constants.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"
#include "util.hpp"

using namespace std::placeholders;


FreddiEvolution::FreddiEvolution(const FreddiArguments &args):
		FreddiState(args, std::bind(&FreddiEvolution::wunction, this, _1, _2, _3, _4)) {}


void FreddiEvolution::step(const double tau) {
	FreddiState::step(tau);
	truncateInnerRadius();
	nonlinear_diffusion_nonuniform_wind_1_2(
			args().calc->tau, args().calc->eps,
			F_in(), Mdot_out(),
			windA(), windB(), windC(),
			wunc(),
			h(), current_.F,
			first(), last());
	truncateOuterRadius();
}


void FreddiEvolution::truncateOuterRadius() {
	if (args().disk->Thot <= 0. ){
		return;
	}
	if (Mdot_in() > Mdot_in_prev()) {
		return;
	}

	auto ii = last() + 1;
	if (args().disk->boundcond == "Teff") {
		do {
			ii--;
		} while( Tph().at(ii) < args().disk->Thot );
	} else if (args().disk->boundcond == "Tirr") {
		do {
			ii--;
		} while( Tirr().at(ii) < args().disk->Thot );
	} else{
		throw std::logic_error("Wrong boundcond");
	}

	if ( ii <= last() - 1 ){
		current_.last = ii;
	}
}


vecd FreddiEvolution::wunction(const vecd &h, const vecd &F, size_t _first, size_t _last) const {
	vecd W(_last + 1, 0.);
	for ( size_t i = _first; i <= _last; ++i ){
		W[i] = pow(std::abs(F[i]), 1. - oprel().m) * pow(h[i], oprel().n) / (1. - oprel().m) / oprel().D;
	}
	return W;
};


FreddiNeutronStarEvolution::NeutronStarStructure::NeutronStarStructure(
		const NeutronStarArguments &args_ns, FreddiEvolution* evolution):

		args_ns(args_ns),
		xi_pow_minus_7_2(std::pow(xi, -3.5)),
		R_x(args_ns.Rx),
		R_m_min(std::max(R_x, evolution->R()[evolution->first()])),
		mu_magn(0.5 * args_ns.Bx * m::pow<3>(R_x)),
		R_dead(args_ns.Rdead > 0. ? args_ns.Rdead : INFINITY),
		F_dead(k_t * xi_pow_minus_7_2 * m::pow<2>(mu_magn) / m::pow<3>(R_dead)),
		R_cor(std::cbrt(evolution->GM() / m::pow<2>(2*M_PI * args_ns.freqx))),
		inverse_beta(args_ns.inversebeta),
		epsilon_Alfven(args_ns.epsilonAlfven),
		hot_spot_area(args_ns.hotspotarea),
		Fmagn(initialize_Fmagn(evolution)),
		dFmagn_dh(initialize_dFmagn_dh(evolution)),
		d2Fmagn_dh2(initialize_d2Fmagn_dh2(evolution)) {
	if (args_ns.Rdead > 0. && args_ns.Rdead < R_cor) {
		throw std::logic_error("R_dead is positive and less than R_cor, it is obvious");
	}
}


vecd FreddiNeutronStarEvolution::NeutronStarStructure::initialize_Fmagn(FreddiEvolution* evolution) const {
	vecd Fmagn_(evolution->Nx());
	const double k = inverse_beta * m::pow<2>(mu_magn) / 3.;
	double brackets;
	for (size_t i = 0; i < evolution->Nx(); i++) {
		if (evolution->R()[i] < R_cor) {
			brackets = -1. + 2. * std::pow(evolution->R()[i] / R_cor, 1.5) - 2./3. * m::pow<3>(evolution->R()[i] / R_cor);
		} else {
			brackets = 1. - 2./3. * std::pow(R_cor / evolution->R()[i], 1.5);
		}
		Fmagn_[i] = k / m::pow<3>(evolution->R()[i]) * brackets;
	}
	return Fmagn_;
}


vecd FreddiNeutronStarEvolution::NeutronStarStructure::initialize_dFmagn_dh(FreddiEvolution* evolution) const {
	vecd dFmagn_dh_(evolution->Nx());
	const double k = inverse_beta * 2 * m::pow<2>(mu_magn) * m::pow<3>(evolution->GM());
	double brackets;
	for (size_t i = 0; i < evolution->Nx(); i++) {
		if (evolution->R()[i] < R_cor) {
			brackets = 1. - std::pow(evolution->R()[i] / R_cor, 1.5);
		} else {
			brackets = -1. + std::pow(R_cor / evolution->R()[i], 1.5);
		}
		dFmagn_dh_[i] = k / m::pow<7>(evolution->h()[i]) * brackets;
	}
	return dFmagn_dh_;
}


vecd FreddiNeutronStarEvolution::NeutronStarStructure::initialize_d2Fmagn_dh2(FreddiEvolution* evolution) const {
	vecd d2Fmagn_dh2_(evolution->Nx());
	const double k = inverse_beta * 2 * m::pow<2>(mu_magn) / evolution->GM();
	double brackets;
	for (size_t i = 0; i < evolution->Nx(); i++) {
		if (evolution->R()[i] < R_cor) {
			brackets = -7. + 4. * std::pow(evolution->R()[i] / R_cor, 1.5);
		} else {
			brackets = 7. - 10. * std::pow(R_cor / evolution->R()[i], 1.5);
		}
		d2Fmagn_dh2_[i] = k / m::pow<4>(evolution->R()[i]) * brackets;
	}
	return d2Fmagn_dh2_;
}


FreddiNeutronStarEvolution::BasicNSMdotFraction::~BasicNSMdotFraction() {}


double FreddiNeutronStarEvolution::NoOutflowNSMdotFraction::operator()(double R_to_Rcor) {
	return 1.;
}


double FreddiNeutronStarEvolution::PropellerNSMdotFraction::operator()(double R_to_Rcor) {
	return 0.;
}


// https://arxiv.org/pdf/1010.1528.pdf Eksi-Kutlu (2010)
// fastness = (R_in/R_cor)^(3/2)  */
double FreddiNeutronStarEvolution::EksiKultu2010NSMdotFraction::operator()(double R_to_Rcor) {
	const double fastness2 = m::pow<3>(R_to_Rcor);
	double p = (1. - 1./fastness2);
	if (p < 0){
		p = 0;
	}
	return 1. - 1.5 * std::sqrt(p) + 0.5 * std::pow(p, 1.5);
}


FreddiNeutronStarEvolution::FreddiNeutronStarEvolution(const FreddiNeutronStarArguments &args):
		FreddiEvolution(args),
		ns_str_(new NeutronStarStructure(*args.ns, this)) {
	// Change initial condition due presence of magnetic field torque. It can spoil user-defined initial disk
	// parameters, such as mass or Fout
	initializeNsMdotFraction();
	if (inverse_beta() <= 0.) {  // F_in is non-zero, Fmang is zero everywhere
		current_.F_in = k_t() * m::pow<2>(mu_magn()) / m::pow<3>(R_cor());
		for (size_t i = 0; i < Nx(); i++) {
			current_.F[i] += current_.F_in;
		}
	} else {  // F_in is zero, and F + Fmagn = initial_cond + Fmang_in
		for (size_t i = 0; i < Nx(); i++) {
			current_.F[i] += -Fmagn()[i] + Fmagn()[0];
		}
	}
}


void FreddiNeutronStarEvolution::initializeNsMdotFraction() {
	if (ns_str_->args_ns.fptype == "no-outflow") {
		fp_.reset(static_cast<BasicNSMdotFraction *>(new NoOutflowNSMdotFraction));
	} else if (ns_str_->args_ns.fptype == "propeller") {
		fp_.reset(static_cast<BasicNSMdotFraction *>(new PropellerNSMdotFraction));
	} else if (ns_str_->args_ns.fptype == "eksi-kultu2010") {
		fp_.reset(static_cast<BasicNSMdotFraction *>(new EksiKultu2010NSMdotFraction));
	} else {
		throw std::logic_error("Wrong fptype");
	}
}


double FreddiNeutronStarEvolution::Lbol_ns() const {
	return eta_ns() * fp() * Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
}


double FreddiNeutronStarEvolution::T_hot_spot() const {
	return std::pow(Lbol_ns() / (4*M_PI * hot_spot_area() * m::pow<2>(R_x()) * GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT), 0.25);
}


double FreddiNeutronStarEvolution::Lx_ns() {
	if (!ns_opt_str_.Lx_ns) {
		const double intensity = Spectrum::Planck_nu1_nu2(T_hot_spot(), args().flux->emin, args().flux->emax, 1e-4);
		ns_opt_str_.Lx_ns = 4*M_PI * hot_spot_area() * m::pow<2>(R_x()) * M_PI * intensity;
	}
	return *ns_opt_str_.Lx_ns;
}


void FreddiNeutronStarEvolution::invalidate_optional_structure() {
	FreddiEvolution::invalidate_optional_structure();
	ns_opt_str_ = NeutronStarOptionalStructure();
}


void FreddiNeutronStarEvolution::truncateInnerRadius() {
	if (R_dead() <= 0.) {
		return;
	}
	if ( Mdot_in() > Mdot_in_prev() ) {
		return;
	}

	const double R_alfven = epsilon_Alfven() *
			std::pow(m::pow<4>(mu_magn()) / (m::pow<2>(Mdot_in()) * GM()), 1./7.);
	double R_m = std::max(R_m_min(), R_alfven);
	R_m = std::min(R_m, R_dead());
	size_t ii;
	for (ii = first(); ii <= last() - 2; ii++) {
		if (R()[ii + 1] > R_m){
			break;
		}
	}
	if (ii >= last() - 2) {
		throw std::runtime_error("R_in > R_out");
	}
	R_m = R()[ii];
	current_.first = ii;

	double new_F_in = 0;
	if (inverse_beta() <= 0.) {
		if (R_m <= R_cor()) {
			new_F_in = k_t() * xi_pow_minus_7_2() * m::pow<2>(mu_magn()) / m::pow<3>(R_cor());
		} else {
			new_F_in = F_dead() * m::pow<3>(R_dead() / R_m);
		}
	}
	current_.F_in = new_F_in;
}


double FreddiNeutronStarEvolution::Mdot_in() const {
	const double dF_dh = (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
	return dF_dh + dFmagn_dh()[first()];
}


vecd FreddiNeutronStarEvolution::windC() const {
	auto C = FreddiEvolution::windC();
	for (size_t i = 0; i < C.size(); i++) {
		C[i] += d2Fmagn_dh2()[i];
	}
	return C;
}


const vecd& FreddiNeutronStarEvolution::Qx() {
	// TODO: do all
	if (!opt_str_.Qx) {
		vecd x(Nx());
		const vecd& CirrCirr = Cirr();
		const double L_disk = (F()[first()] + 0.5 * Mdot_in() * h()[first()]) * omega_i(first());
		const double L_ns = Lbol_ns();
		for (size_t i = first(); i <= last(); i++) {
			x[i] = CirrCirr[i] * (L_ns + L_disk) / (4. * M_PI * m::pow<2>(R()[i]));
		}
		opt_str_.Qx = std::move(x);
	}
	return *opt_str_.Qx;
}


double FreddiNeutronStarEvolution::eta_ns() const {
	const double R_g = GM() / m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
	const double eta = R_g * (1. / R_x() - 0.5 / R()[first()]);
	return eta;
}
