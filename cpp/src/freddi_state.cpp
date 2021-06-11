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
		R_g(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * args.basic->Mx / m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT)),
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
			throw std::invalid_argument("Wrong gridscale");
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
		t(str.args.calc->init_time),
		i_t(0),
		first(0),
		last(str.Nx - 1),
		F(initializeF(str)),
		F_in(0) {}

vecd FreddiState::CurrentState::initializeF(const DiskStructure& str) {
	return str.args.disk->initial_F(str.h);
}


FreddiState::FreddiState(const FreddiArguments& args, const wunc_t& wunc):
		str_(new DiskStructure(args, wunc)),
		current_(*str_),
		disk_irr_source_(initializeFreddiIrradiationSource(args.irr->angular_dist_disk)),
		star_roche_lobe_(str_->semiaxis, args.basic->Mopt / args.basic->Mx, args.basic->roche_lobe_fill),
		star_({}, args.basic->Topt, star_roche_lobe_, args.calc->starlod) {
	initializeWind();
}


FreddiState::FreddiState(const FreddiState& other):
		str_(other.str_),
		current_(other.current_),
		opt_str_(other.opt_str_),
		wind_(other.wind_->clone()),
		disk_irr_source_(other.disk_irr_source_),
		star_roche_lobe_(other.star_roche_lobe_),
		star_(other.star_) {}


void FreddiState::initializeWind() {
	if (args().disk->wind == "no") {
		wind_.reset(static_cast<BasicWind*>(new NoWind(*this)));
	} else if (args().disk->wind == "SS73C") {
		wind_.reset(static_cast<BasicWind*>(new SS73CWind(*this)));
	} else if (args().disk->wind == "Cambier2013") { // Cambier & Smith 1303.6218
		wind_.reset(static_cast<BasicWind*>(new Cambier2013Wind(*this)));
	} else if (args().disk->wind == "__testA__") {
		wind_.reset(static_cast<BasicWind*>(new testAWind(*this)));
	} else if (args().disk->wind == "__testB__") {
		wind_.reset(static_cast<BasicWind*>(new testBWind(*this)));
	} else if (args().disk->wind == "__testC__") {
		wind_.reset(static_cast<BasicWind*>(new testCWind(*this)));
	} else if (args().disk->wind == "__testC_q0_Shields1986__") {
		wind_.reset(static_cast<BasicWind*>(new testCq0Shields1986Wind(*this)));
    } else if (args().disk->wind == "Shields1986") {
        wind_.reset(static_cast<BasicWind*>(new Shields1986Wind(*this)));
    } else if (args().disk->wind == "Janiuk2015") {
        wind_.reset(static_cast<BasicWind*>(new Janiuk2015Wind(*this)));
    } else if (args().disk->wind == "Woods1996AGN") {
        wind_.reset(static_cast<BasicWind*>(new Woods1996Wind(*this)));
    } else if (args().disk->wind == "Woods1996") {
        wind_.reset(static_cast<BasicWind*>(new Woods1996ShieldsApproxWind(*this)));
    } else if (args().disk->wind == "toy") {
        wind_.reset(static_cast<BasicWind*>(new PeriodPaperWind(*this)));
	} else {
		throw std::invalid_argument("Wrong wind");
	}
}


std::shared_ptr<FreddiState::BasicFreddiIrradiationSource> FreddiState::initializeFreddiIrradiationSource(const std::string& angular_dist_type) {
	if (angular_dist_type == "isotropic") {
		return std::make_shared<IsotropicFreddiIrradiationSource>();
	}
	if (angular_dist_type == "plane") {
		return std::make_shared<PlaneFreddiIrradiationSource>();
	}
	throw std::invalid_argument("Wrong angular distribution type");
}


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
}


double FreddiState::Mdot_in() const {
	return (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
}


double FreddiState::Lbol_disk() const {
	return eta() * Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
}


double FreddiState::phase_opt() const {
	return 2.0 * M_PI * (t() - args().flux->ephemeris_t0) / args().basic->period;
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
		x.resize(Nx(), 0.0);
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
		const vecd& K = Kirr();
		const vecd& H = Height();
		const double Lbol = Lbol_disk();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = K[i] * Lbol * angular_dist_disk(H[i] / R()[i]) / (4. * M_PI * m::pow<2>(R()[i]));
		}
		opt_str_.Qx = std::move(x);
	}
	return *opt_str_.Qx;
}


const vecd& FreddiState::Kirr() {
	if(!opt_str_.Kirr) {
		vecd x(Nx());
		const vecd& H = Height();
		for (size_t i = first(); i <= last(); i++) {
			x[i] = args().irr->Cirr * std::pow(H[i] / (R()[i] * 0.05), args().irr->irrindex);
		}
		for (size_t i = last() + 1; i < Nx(); i++) {
			x[i] = args().irr->Cirr_cold * std::pow(H[i] / (R()[i] * 0.05), args().irr->irrindex_cold);
		}
		opt_str_.Kirr = std::move(x);
	}
	return *opt_str_.Kirr;
}


const vecd& FreddiState::Height() {
	if (!opt_str_.Height) {
		vecd x(Nx());
		for (size_t i = first(); i <= last(); i++) {
			x[i] = oprel().Height(R()[i], F()[i]);
		}
		for (size_t i = last() + 1; i < Nx(); i++) {
			x[i] = args().irr->height_to_radius_cold * R()[i];
		}
		opt_str_.Height = std::move(x);
	}
	return *opt_str_.Height;
}


const vecd& FreddiState::Tph_vis() {
	if (!opt_str_.Tph_vis) {
		vecd x(Nx(), 0.0);
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
		vecd x(Nx(), 0.0);
		const double Mdot = std::fabs((F()[first()+1] - F()[first()]) / (h()[first()+1] - h()[first()]));
		for (size_t i = first(); i <= last(); i++) {
			//

			// Qvis due to non-zero Fin:
			x[i] =  3. / (8. * M_PI) * m::pow<4>(GM()) / m::pow<7>(h()[i]) * F()[first()];

			// Qvis due to non-zero Mdot:  = sigma * Trel(dotM)^4
			//      assume that Mdot ~= const where X-rays are generated
			//      dotM = dF/dh
			x[i] += m::pow<4>(Spectrum::T_GR(R()[i], args().basic->kerr, args().basic->Mx, Mdot)) * GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT;

			x[i] = args().flux->colourfactor * std::pow( x[i] / GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT , 0.25);
		}
		opt_str_.Tph_X = std::move(x);
	}
	return *opt_str_.Tph_X;
}


/*
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
*/


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
	IrradiatedStar::sources_t sources;
	sources.push_back(disk_irr_source_->irr_source(*this, Lbol_disk()));
	return sources;
}


double FreddiState::flux_star(const double lambda, const double phase) {
	return star_.luminosity({inclination(), phase}, lambda) / (FOUR_M_PI * m::pow<2>(distance()));
}

double FreddiState::flux_star(const Passband& passband, const double phase) {
	return star_.luminosity({inclination(), phase}, passband) / (FOUR_M_PI * m::pow<2>(distance()));
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
	return lazy_integrate<HotRegion>(opt_str_.Mdot_wind, h(), dMdot_dh);
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
	const double L_edd = 4. * M_PI * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION;;
	const double Mdot_crit = L_edd / (m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta());
	const double eta = 0.025 * 33 * disk->Mdot0 / Mdot_crit;
	const double R_iC = R_IC2out * state.R().back();
	for (size_t i = 0; i < state.Nx(); ++i) {
		const double xi = state.R()[i] / R_iC;
		const double C0 = m_ch0 * (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
		C_[i] = C0 * std::pow((1 + m::pow<2>(0.125 * eta / xi))
								  / (1 + 1. / (m::pow<8>(eta) * m::pow<2>(1 + 262 * m::pow<2>(xi)))), 1. / 6.)
					* std::exp(-m::pow<2>(1 - 1. / std::sqrt(1 + 0.25 / (m::pow<2>(xi)))) / (2 * xi));
	}
}

FreddiState::testAWind::testAWind(const FreddiState& state):
		BasicWind(state),
		kA(state.args().disk->windparams.at("kA")) {
	const double A0 = -kA / m::pow<2>(state.h().back() - state.h().front());
	for (size_t i = 0; i < state.Nx(); ++i) {
		A_[i] = A0 * (state.h()[i] - state.h().front());
	}
}

FreddiState::testBWind::testBWind(const FreddiState& state):
		BasicWind(state),
		kB(state.args().disk->windparams.at("kB")) {
	const double B0 = -kB / m::pow<2>(state.h().back() - state.h().front());
	for (size_t i = 0; i < state.Nx(); ++i) {
		B_[i] = B0;
	}
}

FreddiState::testCWind::testCWind(const FreddiState& state):
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

FreddiState::testCq0Shields1986Wind::testCq0Shields1986Wind(const FreddiState& state):
		BasicWind(state),
		kC(state.args().disk->windparams.at("kC")),
		R_windmin2out(state.args().disk->windparams.at("Rwind")) {}

void FreddiState::testCq0Shields1986Wind::update(const FreddiState& state) {
	BasicWind::update(state);
	const double h_wind_min = std::sqrt(R_windmin2out) * state.h().back();
	for (size_t i = state.first(); i <= state.last(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = -0.5/M_PI * kC * state.Mdot_in() /
					(std::log(1 / R_windmin2out) * m::pow<2>(state.R()[i])) *
					(4 * M_PI * m::pow<3>(state.h()[i])) / m::pow<2>(state.GM());
		}
	}
}

FreddiState::Shields1986Wind::Shields1986Wind(const FreddiState& state):
        BasicWind(state),
	f_X(state.args().disk->windparams.at("f_X")),
	X_f(state.args().disk->windparams.at("X_f")),
        T_iC(state.args().disk->windparams.at("T_iC")) {
    update(state);
}

void FreddiState::Shields1986Wind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    //  1983ApJ...271...70B page 4
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    //  1983ApJ...271...70B page 3
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_iC);
    const double P_0 = (f_X/X_f)*L/(4 * M_PI * GSL_CONST_CGSM_SPEED_OF_LIGHT * m::pow<2>(R_iC));
    const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_iC)/(disk->mu * GSL_CONST_CGSM_MASS_PROTON));
    const double m_ch0 = P_0/C_ch;
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double L_crit = 0.03 * std::sqrt(1e8/T_iC) * L_edd;
    const double el = L/L_crit;

	//std::cerr << "\t" << L << "\t" << R_iC  << "\t" << L_edd << "\t" << C_ch << "\t" << P_0 << "\t" << m_ch0 << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        //  1986ApJ...306...90S page 2
        if (state.R()[i] > 0.1*R_iC) {
            const double xi = state.R()[i] / R_iC;
            const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
            //  1986ApJ...306...90S appendix B page 16
            const double y = std::sqrt(1 + 1/(4*m::pow<2>(xi)) + ((m::pow<2>(xi))/(1 + m::pow<2>(xi)))*((1.2*xi/(xi + el) + 2.2/(1 + m::pow<2>(el)*xi))*(1.2*xi/(xi + el) + 2.2/(1 + m::pow<2>(el)*xi))));
            const double Mach_cc = std::cbrt(((1 + (el + 1)/xi)/(1 + 1/((1 + m::pow<2>(xi))*m::pow<4>(el)))));
	    const double p_po = 0.5 * std::exp(-(((1 - 1/y)*(1 - 1/y))/(2*xi)));
	    C_[i] = -2.0 * C0 * m_ch0 * Mach_cc * p_po * m::pow<2>(y) * std::pow(xi, -5/3);		}
    }
}

FreddiState::Janiuk2015Wind::Janiuk2015Wind(const FreddiState& state):
        BasicWind(state),
        A_0(state.args().disk->windparams.at("A_0")),
        B_1(state.args().disk->windparams.at("B_1")) {
    update(state);
}

void FreddiState::Janiuk2015Wind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double R_g = 2*state.GM()/(m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT));
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double lol = L/L_edd;
	//std::cerr << R_g << "\t" << L/L_edd << "\t" << state.R()[state.last()]/R_g << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        //  https://arxiv.org/pdf/1411.4434.pdf
        if (state.R()[i] > 70.0*R_g) {
			const double fout = 1 - 1/(1+A_0*m::pow<2>(lol));
            const double C_0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
			//C_[i] =  B_1 * 3000.0 * v_r * R_g * fout * ( std::pow((state.R()[state.last()]/R_g), 0.2) - std::pow((state.R()[i]/R_g), 0.2) );
            //*( std::pow((state.R()[state.last()]/R_g), 0.2) - std::pow((state.R()[i]/R_g), 0.2)
            const double Q = 2 * (3/(8 * M_PI)) * ((m::pow<4>(state.GM()))/(m::pow<7>(state.h()[i])))  ;
            B_[i] =  - 0.75 * C_0 * (1/B_1) * ((4 * Q * state.R()[i])/(3* state.GM()))*fout ; }
    }
}

FreddiState::Woods1996Wind::Woods1996Wind(const FreddiState& state):
        BasicWind(state),
        C_0(state.args().disk->windparams.at("C_0")),
        T_iC(state.args().disk->windparams.at("T_iC")) {
    update(state);
}

void FreddiState::Woods1996Wind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_iC);
    const double le = L/L_edd;
    double R_tr;
    if ( le <= 0.01 ) {
        R_tr = 6.0*R_iC;
    }
    else
        {
            R_tr = R_iC*(6.0 + 5.4*std::log(le/0.01) + 4.1*((std::log(le/0.01))*(std::log(le/0.01))));
    }

    double f_L;
    if ( le <= 0.1 ) {
        f_L = 1.0 ;
    }
    else {
            f_L = std::pow((0.1/le), 0.15) ;
    }


    //std::cerr << "\t" << L << "\t" << R_iC  << "\t" << L_edd << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        //
        double g_R;
        if ( state.R()[i] <= R_tr ) {
            g_R = 1.0 ;
        }
        else {
            g_R = state.R()[i]/R_tr ;
        }
        const double xi1 = R_iC / state.R()[i];
        const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
        const double ExP = std::exp(-(((1.0 - (1/std::sqrt(1.0 + 0.25*m::pow<2>(xi1))))*(1.0 - (1/std::sqrt(1.0 + 0.25*m::pow<2>(xi1)))))/(2.0/xi1)));
        C_[i] = -2.0 *(2e42/(state.args().basic->Mx))*(C_0/1e13) * C0 * le * m::pow<2>(xi1) * f_L * g_R * ExP ;		}
    }



FreddiState::Woods1996ShieldsApproxWind::Woods1996ShieldsApproxWind(const FreddiState& state):
BasicWind(state),
        Xi_max(state.args().disk->windparams.at("Xi_max")),
        T_iC(state.args().disk->windparams.at("T_iC")),
        W_pow(state.args().disk->windparams.at("W_pow")) {
    update(state);
}

void FreddiState::Woods1996ShieldsApproxWind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_iC);
    const double VeL = std::sqrt(state.GM()/R_iC) ;
    const double C_iC = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_iC)/( GSL_CONST_CGSM_MASS_PROTON));
    //const double m_ch0 = disk->Mdot0 / (M_PI * state.R().back()*state.R().back());
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double L_crit = 0.03 * std::sqrt(1e8/T_iC) * L_edd;
    const double el = L / L_crit;



    //std::cerr << "\t" << L << " \t" << L_edd << "\t" << L/L_edd  << "\t" << L_crit << "\t" << R_iC << "\t" << state.eta() << "\t" << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        if (state.R()[i] > 0.1*R_iC) {
            //  1986ApJ...306...90S page 2
            const double xi = state.R()[i] / R_iC;
            const double xi1 = R_iC / state.R()[i];
            const double T_ch = T_iC * std::pow(el, 2.0 / 3.0) * std::pow(xi, -2.0 / 3.0);
            const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ch) / (disk->mu * GSL_CONST_CGSM_MASS_PROTON));
            const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
            const double Fr =
                    L / (4.0 * M_PI * m::pow<2>(state.R()[i]) * Xi_max * C_ch * GSL_CONST_CGSM_SPEED_OF_LIGHT);
//const double Fc = std::pow(((1.0 + ( ((0.125 * el + 0.00382)/ xi) *((0.125 * el + 0.00382)/ xi) ))/( 1 + 1/( (el*el*el*el*(1 + 262.0*xi*xi))*(el*el*el*el*(1 + 262.0*xi*xi)) ) ) ), 1.0/6.0) ;
            const double Fc = ((std::pow((1 + m::pow<2>(((0.125 * el + 0.00382) * xi1))), (1.0 / 6.0))) /
                               std::pow((1 + m::pow<-2>((m::pow<4>(el)* (1.0 + 262.0 * m::pow<2>(xi))))),
                                        (1.0 / 6.0)));
	    //const double Fc = ((std::pow((1 + std::pow(((0.125 * el + 0.00382) * xi1), 2.0)), (1.0 / 6.0))) /
              //                 std::pow((1 + std::pow((el * el * el * el * (1.0 + 262.0 * xi * xi)), -2.0)),
                //                        (1.0 / 6.0)));
            const double Expo = std::exp(-(((1.0 - (1 / std::sqrt(1.0 + 0.25 * m::pow<2>(xi1)))) *
                                            (1.0 - (1 / std::sqrt(1.0 + 0.25 * m::pow<2>(xi1))))) / (2.0 * xi)));
            C_[i] = - 2.0 * W_pow * C0 * Fr * Fc * Expo;
        }

    }
}

FreddiState::PeriodPaperWind::PeriodPaperWind(const FreddiState& state):
	BasicWind(state),
	W_pow(state.args().disk->windparams.at("W_pow")) {
    update(state);
}

void FreddiState::PeriodPaperWind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;

    for (size_t i = state.first(); i <= state.last(); ++i) {
	//const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
	const double Mdot = state.Mdot_in() * (state.h()[i] - state.h()[state.first()]) /(m::pow<2>(state.h()[state.last()] - state.h()[state.first()])) ;
	C_[i] = - 2.0 * W_pow * Mdot;
    }
}

FreddiState::BasicFreddiIrradiationSource::~BasicFreddiIrradiationSource() {}


Vec3 FreddiState::BasicFreddiIrradiationSource::position(const FreddiState& state) const {
	return {-state.semiaxis(), 0.0, 0.0};
}

double FreddiState::BasicFreddiIrradiationSource::Height2R(FreddiState& state) const {
	const auto& H = state.Height();
	const auto& R = state.R();
	double max_H2R = 0.0;
	for (size_t i = state.first(); i < state.Nx(); i++) {
		max_H2R = std::max(H[i] / R[i], max_H2R);
	}
	return max_H2R;
}

double FreddiState::IsotropicFreddiIrradiationSource::angular_dist(const double mu) const {
	return 1.0;
}

std::unique_ptr<IrrSource> FreddiState::IsotropicFreddiIrradiationSource::irr_source(FreddiState& state, const double luminosity) const {
	return std::make_unique<PointAccretorSource>(position(state), luminosity, state.args().flux->star_albedo, Height2R(state));
}

FreddiState::PlaneFreddiIrradiationSource::PlaneFreddiIrradiationSource():
		normal(0.0, 0.0) {}

double FreddiState::PlaneFreddiIrradiationSource::angular_dist(const double mu) const {
	return 2.0 * mu;
}

std::unique_ptr<IrrSource> FreddiState::PlaneFreddiIrradiationSource::irr_source(FreddiState& state, const double luminosity) const {
	return std::make_unique<CentralDiskSource>(position(state), normal, luminosity, state.args().flux->star_albedo, Height2R(state));
}

double FreddiState::Sigma_minus(double r) const {
	// Lasota et al., A&A 486, 523–528 (2008), Eq A.1, DOI: 10.1051/0004-6361:200809658
	return 74.6 * std::pow(args().basic->alphacold / 0.1, -0.83) * std::pow(r / 1e10, 1.18)
		* std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.40);
}

double FreddiState::Sigma_plus(double r) const {
	// Lasota et al., A&A 486, 523–528 (2008), Eq A.1, DOI: 10.1051/0004-6361:200809658
	return 39.9 * std::pow(args().basic->alpha / 0.1, -0.80) * std::pow(r / 1e10, 1.11)
		* std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.37);
}

double FreddiState::v_cooling_front(double r) {
        // The cooling-front velocity depends on the ratio between the current Sigma and critical Sigmas
        // Ludwig et al., A&A 290, 473-486 (1994), section 3
        // units: cm/s
        const double Sigma_plus_ = Sigma_plus(r);
        const double sigma =  std::log( Sigma()[last()] / Sigma_plus_ ) /  std::log( Sigma_minus(r)/Sigma_plus_ ) ;
        return 1e5 * (1.439-5.305*sigma+10.440*m::pow<2>(sigma)-10.55*m::pow<3>(sigma)+4.142*m::pow<4>(sigma))
               * std::pow(args().basic->alpha / 0.2, 0.85-0.69*sigma) 
               * std::pow(args().basic->alphacold / 0.05, 0.05+0.69*sigma)
               * std::pow(r / 1e10, 0.035)
               * std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.012);
}

double FreddiState::R_cooling_front(double r)  {
        // previous location of Rhot moves with the cooling-front velocity:
        return  R()[last()] - v_cooling_front(r) * args().calc->tau;       
        //return  R()[last()] - v_cooling_front(R()[last()]) * args().calc->tau  ; 
        // this variant leads to more abrupt evolution, since the front velocity is larger
}


