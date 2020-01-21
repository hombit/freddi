#include "freddi_state.hpp"

#include <cmath>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"


FreddiState::DiskStructure::DiskStructure(const FreddiArguments &args, const wunc_t& wunc):
		args(args),
		Nt(static_cast<size_t>(std::round(args.calc->time / args.calc->tau))),
		Nx(args.calc->Nx),
		GM(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args.basic->Mx),
		eta(efficiency_of_accretion(args.basic->kerr)),
		semiaxis(args.basic->semiaxis(args.basic->Mx, args.basic->Mopt, args.basic->period)),
		inclination(args.flux->inclination / 180.0 * M_PI),
		cosi(std::cos(args.flux->inclination / 180.0 * M_PI)),
		distance(args.flux->distance),
		cosiOverD2(cosi / m::pow<2>(distance)),
		oprel(*(args.disk->oprel.get())),
		h(initialize_h(args, Nx)),
		R(initialize_R(h, GM)),
		wunc(wunc) {}

vecd FreddiState::DiskStructure::initialize_h(const FreddiArguments& args, size_t Nx) {
	const double h_in = args.basic->h(args.basic->rin);
	const double h_out = args.basic->h(args.basic->rout);

	vecd h(Nx);
	for (size_t i = 0; i < Nx; i++) {
		if (args.calc->gridscale == "log") {
			h[i] = h_in * pow(h_out / h_in, i / (Nx - 1.));
		} else if (args.calc->gridscale == "linear") {
			h[i] = h_in + (h_out - h_in) * i / (Nx - 1.);
		} else {
			throw std::logic_error("Wrong gridscale");
		}
	}
	return h;
}

vecd FreddiState::DiskStructure::initialize_R(const vecd& h, double GM) {
	vecd R(h.size());
	for (size_t i = 0; i < h.size(); i++) {
		R[i] = m::pow<2>(h[i]) / GM;
	}
	return R;
}


FreddiState::CurrentState::CurrentState(const DiskStructure& str):
		Mdot_out(str.args.disk->Mdotout),
		last(str.Nx - 1),
		F(initializeF(str)) {}

vecd FreddiState::CurrentState::initializeF(const DiskStructure& str) {
	const auto& disk = str.args.disk;
	const auto& h = str.h;
	const auto& oprel = str.oprel;

	const double h_in = h.front();
	const double h_out = h.back();

	vecd F(str.Nx);

	if (disk->initialcond == "power" or disk->initialcond == "powerF") {
		for (size_t i = 0; i < str.Nx; ++i) {
			F[i] = disk->F0 * std::pow((h[i] - h_in) / (h_out - h_in), disk->powerorder);
		}
	} else if (disk->initialcond == "powerSigma") {
		for (size_t i = 0; i < str.Nx; ++i) {
			const double Sigma_to_Sigmaout = std::pow((h[i] - h_in) / (h_out - h_in), disk->powerorder);
			F[i] = disk->F0 * std::pow(h[i] / h_out, (3. - oprel.n) / (1. - oprel.m)) *
					std::pow(Sigma_to_Sigmaout, 1. / (1. - oprel.m));
		}
	} else if (disk->initialcond == "sinusF" || disk->initialcond == "sinus") {
		for (size_t i = 0; i < str.Nx; ++i) {
			F[i] = disk->F0 * std::sin((h[i] - h_in) / (h_out - h_in) * M_PI_2);
		}
	} else if (disk->initialcond == "quasistat") {
		for (size_t i = 0; i < str.Nx; ++i) {
			const double xi_LS2000 = h[i] / h_out;
			F[i] = disk->F0 * oprel.f_F(xi_LS2000) * (1. - h_in / h[i]) / (1. - h_in / h_out);
		}
	} else if (disk->initialcond == "gaussF") {
		for (size_t i = 0; i < str.Nx; ++i) {
			const double xi = (h[i] - h_in) / (h_out - h_in);
			F[i] = disk->F0 *
					std::exp(-m::pow<2>(xi - disk->gaussmu) / (2. * m::pow<2>(disk->gausssigma)));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}

	return F;
}


FreddiState::FreddiState(const FreddiArguments& args, const wunc_t& wunc):
		 str_(new DiskStructure(args, wunc)),
		 current_(*str_),
		 star_({}, args.basic->Topt, args.basic->Ropt, args.calc->starlod) {
	initializeWind();
}


void FreddiState::initializeWind() {
	if (args().disk->wind == "no") {
		wind_.reset(static_cast<BasicWind*>(new NoWind(*this)));
	} else if (args().disk->wind == "SS73C") {
		wind_.reset(static_cast<BasicWind*>(new SS73CWind(*this)));
	} else if (args().disk->wind == "Cambier2013") { // Cambier & Smith 1303.6218
		wind_.reset(static_cast<BasicWind*>(new Cambier2013Wind(*this)));
	} else if (args().disk->wind == "__testA__") {
		wind_.reset(static_cast<BasicWind*>(new __testA__Wind(*this)));
	} else if (args().disk->wind == "__testB__") {
		wind_.reset(static_cast<BasicWind*>(new __testB__Wind(*this)));
	} else if (args().disk->wind == "__testC__") {
		wind_.reset(static_cast<BasicWind*>(new __testC__Wind(*this)));
	} else if (args().disk->wind == "__testC_q0_Shields1986__") {
		wind_.reset(static_cast<BasicWind*>(new __testC_q0_Shields1986__(*this)));
	} else {
		throw std::logic_error("Wrong wind");
	}
}


FreddiState::FreddiState(const FreddiState& other):
		str_(other.str_),
		current_(other.current_),
		opt_str_(other.opt_str_),
		wind_(other.wind_->clone()),
		star_(other.star_) {}


void FreddiState::invalidate_optional_structure() {
	opt_str_ = DiskOptionalStructure();
}


void FreddiState::replaceArgs(const FreddiArguments &args) {
	str_.reset(new DiskStructure(args, wunc()));
	invalidate_optional_structure();
}


void FreddiState::step(double tau) {
	set_Mdot_in_prev();
	invalidate_optional_structure();
	current_.i_t ++;
	current_.t += tau;
	wind_->update(*this);
	star_.set_sources(star_irr_sources());
}


double FreddiState::Mdot_in() const {
	return (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
}


double FreddiState::Lbol_disk() const {
	return eta() * Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
}


double FreddiState::Lx() {
	if (!opt_str_.Lx) {
		opt_str_.Lx = Luminosity(Tph_X(), args().flux->emin, args().flux->emax) / m::pow<4>(args().flux->colourfactor);
	}
	return *opt_str_.Lx;
}


const vecd& FreddiState::W() {
	if (!opt_str_.W) {
		auto x = wunc()(h(), F(), first(), last());
		x.resize(Nx());
		opt_str_.W = std::move(x);
	}
	return *opt_str_.W;
}


const vecd& FreddiState::Sigma() {
	if (!opt_str_.Sigma) {
		vecd x(Nx());
		const vecd& WW = W();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = WW[i] * m::pow<2>(GM()) / (4. * M_PI * m::pow<3>(h()[i]));
		}
		opt_str_.Sigma = std::move(x);
	}
	return *opt_str_.Sigma;
}


const vecd& FreddiState::Tph() {
	if (!opt_str_.Tph) {
		vecd x(Nx());
		const vecd& Tvis = Tph_vis();
		const vecd& QxQx = Qx();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = std::pow(m::pow<4>(Tvis[i]) + QxQx[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		opt_str_.Tph = std::move(x);
	}
	return *opt_str_.Tph;
}


const vecd& FreddiState::Tirr() {
	if (!opt_str_.Tirr) {
		vecd x(Nx());
		const vecd& QxQx = Qx();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = std::pow(QxQx[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		opt_str_.Tirr = std::move(x);
	}
	return *opt_str_.Tirr;
}


const vecd& FreddiState::Qx() {
	if (!opt_str_.Qx) {
		vecd x(Nx());
		const vecd& CirrCirr = Cirr();
		const double Lbol = Lbol_disk();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = CirrCirr[i] * Lbol / (4. * M_PI * m::pow<2>(R()[i]));
		}
		opt_str_.Qx = std::move(x);
	}
	return *opt_str_.Qx;
}


const vecd& FreddiState::Cirr() {
	if (!opt_str_.Cirr) {
		if (args().irr->irrfactortype == "const") {
			opt_str_.Cirr = vecd(Nx(), args().irr->Cirr);
		} else if (args().irr->irrfactortype == "square") {
			vecd x(Nx());
			const vecd& H = Height();
			for (size_t i = first(); i < Nx(); i++) {
				x[i] = args().irr->Cirr * (H[i] / R()[i]) * (H[i] / R()[i]);
			}
			opt_str_.Cirr = std::move(x);
		} else {
			throw std::logic_error("Wrong irrfactor");
		}
	}
	return *opt_str_.Cirr;
}


const vecd& FreddiState::Height() {
	if (!opt_str_.Height) {
		vecd x(Nx());
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = oprel().Height(R()[i], F()[i]);
		}
		opt_str_.Height = std::move(x);
	}
	return *opt_str_.Height;
}


const vecd& FreddiState::Tph_vis() {
	if (!opt_str_.Tph_vis) {
		vecd x(Nx());
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = (GM() * std::pow(h()[i], -1.75)
					* std::pow(3. / (8. * M_PI) * F()[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25));
		}
		opt_str_.Tph_vis = std::move(x);
	}
	return *opt_str_.Tph_vis;
}


const vecd& FreddiState::Tph_X() {
	if (!opt_str_.Tph_X) {
		vecd x(Nx());
		x[first()] = 0;
		for (size_t i = first() + 1; i < Nx(); i++) {
			x[i] = (args().flux->colourfactor * Spectrum::T_GR(R()[i], args().basic->kerr, args().basic->Mx,
					F()[i] / (h()[i] - h()[first()])));
		}
		opt_str_.Tph_X = std::move(x);
	}
	return *opt_str_.Tph_X;
}


double FreddiState::lazy_magnitude(boost::optional<double>& m, double lambda, double F0) {
	if (!m) {
		m = magnitude(lambda, F0);
	}
	return *m;
}


double FreddiState::Luminosity(const vecd& T, double nu1, double nu2) const {
	// 2 - two sides
	// pi = \int cos(phi) dtheta dphi
	// pi * \int Bnu dnu gives flux
	// 2 * \int dS gives disk area on both sides
	return 2. * M_PI * integrate<HotRegion>([&T, nu1, nu2](const size_t i) -> double { return Spectrum::Planck_nu1_nu2(T[i], nu1, nu2, 1e-4); });
}


IrradiatedStar::sources_t FreddiState::star_irr_sources() {
	const Vec3 position(-semiaxis(), 0.0, 0.0);
	const UnitVec3 normal(0.0, 0.0);
	const double luminosity = Lbol_disk();
	const double Height2R = Height()[last()] / R()[last()];

	IrradiatedStar::sources_t sources;
	sources.push_back(std::make_unique<CentralDiskSource>(position, normal, luminosity, Height2R));
	return sources;
}


double FreddiState::flux_star(const double lambda, const double phase) {
	return star_.dLdOmega({inclination(), phase}, lambda) / m::pow<2>(distance());
}

double FreddiState::flux_star(const Passband& passband, const double phase) {
	return star_.dLdOmega({inclination(), phase}, passband) / m::pow<2>(distance());
}


FreddiState::BasicWind::BasicWind(const FreddiState &state):
		A_(state.Nx(), 0.), B_(state.Nx(), 0.), C_(state.Nx(), 0.) {}

FreddiState::BasicWind::~BasicWind() = default;

FreddiState::SS73CWind::SS73CWind(const FreddiState &state):
		BasicWind(state) {
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state.GM();
	const double Mdot_crit = L_edd / (m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta());
	for (size_t i = 0; i < state.Nx(); ++i) {
		C_[i] = -Mdot_crit / (2 * M_PI * state.R().front() * state.R()[i]) *
				(4 * M_PI * m::pow<3>(state.h()[i])) / m::pow<2>(state.GM());
	}
}

FreddiState::Cambier2013Wind::Cambier2013Wind(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at("kC")),
		R_IC2out(state.args().disk->windparams.at("RIC")) {
	const auto disk = state.args().disk;
	const double m_ch0 = -kC * disk->Mdot0 / (M_PI * m::pow<2>(state.R().back()));  // dM / dA
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state.GM();
	const double Mdot_crit = L_edd / (m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta());
	const double eta = 0.025 * 33 * disk->Mdot0 / Mdot_crit;
	const double R_iC = R_IC2out * state.R().back();
	for (size_t i = 0; i < state.Nx(); ++i) {
		const double xi = state.R()[i] / R_iC;
		const double C0 = m_ch0 * (4 * M_PI * m::pow<3>(state.h()[i])) / m::pow<2>(state.GM());
		C_[i] = C0 * std::pow((1 + m::pow<2>(0.125 * eta / xi))
								  / (1 + 1. / (m::pow<8>(eta) * m::pow<2>(1 + 262 * m::pow<2>(xi)))), 1. / 6.)
					* std::exp(-m::pow<2>(1 - 1. / std::sqrt(1 + 0.25 / m::pow<2>(xi))) / (2 * xi));
	}
}

FreddiState::__testA__Wind::__testA__Wind(const FreddiState& state):
		BasicWind(state),
		kA(state.args().disk->windparams.at("kA")) {
	const double A0 = -kA / m::pow<2>(state.h().back() - state.h().front());
	for (size_t i = 0; i < state.Nx(); ++i) {
		A_[i] = A0 * (state.h()[i] - state.h().front());
	}
}

FreddiState::__testB__Wind::__testB__Wind(const FreddiState& state):
		BasicWind(state),
		kB(state.args().disk->windparams.at("kB")) {
	const double B0 = -kB / m::pow<2>(state.h().back() - state.h().front());
	for (size_t i = 0; i < state.Nx(); ++i) {
		B_[i] = B0;
	}
}

FreddiState::__testC__Wind::__testC__Wind(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at("kC")) {
	const double C0 = kC * state.args().disk->Mdotout / (state.h().back() - state.h().front());
	const double h_wind_min = state.h().back() / 2;
	for (size_t i = 0; i < state.Nx(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = C0 * 0.5 * (std::cos(2. * M_PI * (state.h()[i] - h_wind_min) / (state.h().back() - h_wind_min)) - 1);
		}
	}
}

FreddiState::__testC_q0_Shields1986__::__testC_q0_Shields1986__(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at("kC")),
		R_windmin2out(state.args().disk->windparams.at("Rwind")) {}

void FreddiState::__testC_q0_Shields1986__::update(const FreddiState& state) {
	BasicWind::update(state);
	const double h_wind_min = std::sqrt(R_windmin2out) * state.h().back();
	for (size_t i = state.first(); i <= state.last(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = -0.5/M_PI * kC * state.Mdot_in() /
					(std::log(1 / R_windmin2out) * state.R()[i] * state.R()[i]) *
					(4 * M_PI * m::pow<3>(state.h()[i])) / m::pow<2>(state.GM());
		}
	}
}
