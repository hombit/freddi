#include "freddi_state.hpp"

#include <cmath>
#include <exception>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "constants.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"


FreddiState::DiskStructure::DiskStructure(const FreddiArguments &args, const wunc_t& wunc):
		args(args),
		Nt(static_cast<size_t>(std::round(args.calc->time / args.calc->tau))),
		Nx(args.calc->Nx),
		GM(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args.basic->Mx),
		eta(efficiency_of_accretion(args.basic->kerr)),
		cosi(std::cos(args.flux->inclination / 180 * M_PI)),
		distance(args.flux->distance),
		cosiOverD2(cosi / (distance*distance)),
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
		R[i] = h[i]*h[i] / GM;
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
					std::exp(-(xi - disk->gaussmu)*(xi - disk->gaussmu) / (2. * disk->gausssigma*disk->gausssigma));
		}
	} else {
		throw std::logic_error("Wrong initialcond");
	}

	return F;
}


FreddiState::FreddiState(const FreddiArguments& args, const wunc_t& wunc):
		 str_(new DiskStructure(args, wunc)),
		 current_(*str_) {
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
    } else if (args().disk->wind == "__Unstedy_Test_Hunter__") {
        wind_.reset(static_cast<BasicWind*>(new __Unstedy_Test_Hunter__(*this)));
	} else {
		throw std::logic_error("Wrong wind");
	}
}


FreddiState::FreddiState(const FreddiState& other):
		str_(other.str_),
		current_(other.current_),
		opt_str_(other.opt_str_),
		wind_(other.wind_->clone()) {}


void FreddiState::step(double tau) {
	set_Mdot_in_prev();
	invalidate_disk_optional_structure();
	current_.i_t ++;
	current_.t += tau;
	wind_->update(*this);
}


double FreddiState::Mdot_in() const {
	return (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
}

double FreddiState::Lx() {
	if (!opt_str_.Lx) {
		opt_str_.Lx = Luminosity(Tph_X(), args().flux->emin, args().flux->emax) * pow(args().flux->colourfactor, -4);
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
		for (size_t i = first(); i <= last(); i++) {
			x[i] = WW[i] * GM() * GM() / (4. * M_PI * pow(h()[i], 3.));
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
		for (size_t i = first(); i <= last(); i++) {
			x[i] = std::pow(std::pow(Tvis[i], 4) + QxQx[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT, 0.25);
		}
		opt_str_.Tph = std::move(x);
	}
	return *opt_str_.Tph;
}


const vecd& FreddiState::Tirr() {
	if (!opt_str_.Tirr) {
		vecd x(Nx());
		const vecd& QxQx = Qx();
		for (size_t i = first(); i <= last(); i++) {
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
		const double Mdot = Mdot_in();
		for (size_t i = first(); i <= last(); i++) {
			x[i] = (CirrCirr[i] * eta() * Mdot * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT
					/ (4. * M_PI * R()[i] * R()[i]));
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
			for (size_t i = first(); i <= last(); i++) {
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
		for (size_t i = first(); i <= last(); i++) {
			x[i] = oprel().Height(R()[i], F()[i]);
		}
		opt_str_.Height = std::move(x);
	}
	return *opt_str_.Height;
}


const vecd& FreddiState::Tph_vis() {
	if (!opt_str_.Tph_vis) {
		vecd x(Nx());
		for (size_t i = first(); i <= last(); i++) {
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
		for (size_t i = first() + 1; i <= last(); i++) {
			x[i] = (args().flux->colourfactor * Spectrum::T_GR(R()[i], args().basic->kerr, args().basic->Mx,
					F()[i] / (h()[i] - h()[first()]), R()[first()]));
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


double FreddiState::I_lambda(const double lambda) {
	const vecd& T = Tph();
	return integrate([&T, lambda](const size_t i) -> double { return Spectrum::Planck_lambda(T[i], lambda); });
}


double FreddiState::Luminosity(const vecd& T, double nu1, double nu2) const {
	return 2. * M_PI * integrate(
			[&T, nu1, nu2](const size_t i) -> double { return Spectrum::Planck_nu1_nu2(T[i], nu1, nu2, 1e-4); }
	);
}


double FreddiState::Mdot_wind() {
	auto dMdot_dh = [this](const size_t i) -> double {
		double dFdh;
		if (i == first()) {
			dFdh = (F()[i+1] - F()[i]) / (h()[i+1] - h()[i]);
		} else if (i == last()) {
			dFdh = (F()[i] - F()[i-1]) / (h()[i] - h()[i-1]);
		} else {
			const double delta_0 = h()[i] - h()[i-1];
			const double delta_1 = h()[i+1] - h()[i];
			dFdh = (F()[i+1] * delta_0 * delta_0 / (delta_0 + delta_1) +
					F()[i] * (delta_1 - delta_0) -
					F()[i-1] * delta_1 * delta_1 / (delta_0 + delta_1)) /
					(delta_0 * delta_1);
		}
		// Wind loss rate sign is opposite disk loss rate sign, e.g. usually it should be positive
		return -(windA()[i] * dFdh + windB()[i] * F()[i] + windC()[i]);
	};
	return lazy_integrate(opt_str_.Mdot_wind, h(), dMdot_dh);
}



FreddiState::BasicWind::BasicWind(const FreddiState &state):
		A_(state.Nx(), 0.), B_(state.Nx(), 0.), C_(state.Nx(), 0.) {}

FreddiState::BasicWind::~BasicWind() = default;

FreddiState::SS73CWind::SS73CWind(const FreddiState &state):
		BasicWind(state) {
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state.GM();
	const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT *
			state.eta());
	for (size_t i = 0; i < state.Nx(); ++i) {
		C_[i] = -Mdot_crit / (2 * M_PI * state.R().front() * state.R()[i]) *
				(4 * M_PI * state.h()[i] * state.h()[i] * state.h()[i]) / (state.GM() * state.GM());
	}
}

FreddiState::Cambier2013Wind::Cambier2013Wind(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at(0)),
		R_IC2out(state.args().disk->windparams.at(1)) {
	const auto disk = state.args().disk;
	const double m_ch0 = -kC * disk->Mdot0 / (M_PI * state.R().back()*state.R().back());  // dM / dA
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
						 GSL_CONST_CGSM_THOMSON_CROSS_SECTION * state.GM();
	const double Mdot_crit = L_edd / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * state.eta());
	const double eta = 0.025 * 33 * disk->Mdot0 / Mdot_crit;
	const double R_iC = R_IC2out * state.R().back();
	for (size_t i = 0; i < state.Nx(); ++i) {
		const double xi = state.R()[i] / R_iC;
		const double C0 = m_ch0 * (4 * M_PI * state.h()[i]*state.h()[i]*state.h()[i]) / (state.GM()*state.GM());
		C_[i] = C0 * std::pow((1 + std::pow(0.125 * eta / xi, 2))
								  / (1 + 1. / (std::pow(eta, 8) * std::pow(1 + 262 * xi * xi, 2))), 1. / 6.)
					* std::exp(-std::pow(1 - 1. / std::sqrt(1 + 0.25 / (xi * xi)), 2) / (2 * xi));
	}
}

FreddiState::__testA__Wind::__testA__Wind(const FreddiState& state):
		BasicWind(state),
		kA(state.args().disk->windparams.at(0)) {
	const double A0 = -kA / ((state.h().back() - state.h().front()) * (state.h().back() - state.h().front()));
	for (size_t i = 0; i < state.Nx(); ++i) {
		A_[i] = A0 * (state.h()[i] - state.h().front());
	}
}

FreddiState::__testB__Wind::__testB__Wind(const FreddiState& state):
		BasicWind(state),
		kB(state.args().disk->windparams.at(0)) {
	const double B0 = -kB / ((state.h().back() - state.h().front()) * (state.h().back() - state.h().front()));
	for (size_t i = 0; i < state.Nx(); ++i) {
		B_[i] = B0;
	}
}

FreddiState::__testC__Wind::__testC__Wind(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at(0)) {
	const double C0 = kC * state.args().disk->Mdotout / (state.h().back() - state.h().front());
	const double h_wind_min = state.h().back() / 2;
	for (size_t i = 0; i < state.Nx(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = C0 * 0.5 *
					(std::cos(2. * M_PI * (state.h()[i] - h_wind_min) / (state.h().back() - h_wind_min)) - 1);
		}
	}
}

FreddiState::__testC_q0_Shields1986__::__testC_q0_Shields1986__(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at(0)),
		R_windmin2out(state.args().disk->windparams.at(1)) {}

void FreddiState::__testC_q0_Shields1986__::update(const FreddiState& state) {
	BasicWind::update(state);
	const double h_wind_min = std::sqrt(R_windmin2out) * state.h().back();
	for (size_t i = state.first(); i <= state.last(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = -0.5/M_PI * kC * state.Mdot_in() /
					(std::log(1 / R_windmin2out) * state.R()[i] * state.R()[i]) *
					(4 * M_PI * state.h()[i]*state.h()[i]*state.h()[i]) / (state.GM()*state.GM());
		}
	}
}
FreddiState::__Unstedy_Test_Hunter__::__Unstedy_Test_Hunter__(const FreddiState& state):
        BasicWind(state),
		fXI(state.args().disk->windparams.at(0)),
        T_iC(state.args().disk->windparams.at(1)) {
	update(state);
}

void FreddiState::__Unstedy_Test_Hunter__::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    //  1983ApJ...271...70B page 4
    const double L = state.Mdot_in() * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT * state.eta();
    //  1983ApJ...271...70B page 3
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_iC);
	const double P_0 = fXI*L/(4 * M_PI * GSL_CONST_CGSM_SPEED_OF_LIGHT * R_iC * R_iC);
	const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_iC)/(disk->mu * GSL_CONST_CGSM_MASS_PROTON));
	const double m_ch0 = P_0/C_ch;
    //const double m_ch0 = disk->Mdot0 / (M_PI * state.R().back()*state.R().back());
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double L_crit = 0.03 * std::sqrt(1e8/T_iC) * L_edd;
	const double el = L/L_crit;

	std::cerr << L << "\t" << L_crit << "\t" << T_iC << "\t" << R_iC  << "\t" << L_edd << "\t" << C_ch << "\t" << m_ch0 << std::endl;

	for (size_t i = state.first(); i <= state.last(); ++i) {
        //  1986ApJ...306...90S page 2
        if (state.R()[i] > 0.1*R_iC) {
            const double xi = state.R()[i] / R_iC;
            const double C0 = (4 * M_PI * state.h()[i]*state.h()[i]*state.h()[i]) / (state.GM()*state.GM());
            //  1986ApJ...306...90S appendix B page 16
            const double y = std::sqrt(1 + 1/(4*xi*xi) + ((xi*xi)/(1 + xi*xi))*((1.2*xi/(xi + el) + 2.2/(1 + el*el*xi))*(1.2*xi/(xi + el) + 2.2/(1 + el*el*xi))));
            const double Mach_cc = std::pow(((1 + (el + 1)/xi)/(1 + 1/((1 + xi*xi)*el*el*el*el))) , 1/3);
			const double p_po = 0.5 * std::exp(-(((1 - 1/y)*(1 - 1/y))/(2*xi)));
			C_[i] = -2.0 * C0 * m_ch0 * Mach_cc * p_po * y * y * std::pow(xi, -5/3);		}
    }
}
