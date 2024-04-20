#include "freddi_state.hpp"

#include <cmath>
#include <string>

#include "gsl_const_cgsm.h"

#include "arguments.hpp"
#include "nonlinear_diffusion.hpp"
#include "orbit.hpp"
#define VERB_LEVEL_MESSAGES 30 

// NOTES: 2 lines marked by // @XRPCALCApril24
//  have to be edited uncommented shadow=0, Tirr_crit=1e10 for scatter=no
// to revive the version used for XRPs in April 2024

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
		oprel(args.disk->oprel),
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
		first(initializeFirst(str)),
		last(str.Nx - 1),
		F(initializeF(str)),
		F_in(0) {}

size_t FreddiState::CurrentState::initializeFirst(const DiskStructure& str) {
	return str.args.disk->initial_first(str.h);
}

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
	} else if (args().disk->wind == "ShieldsOscil1986") {
		wind_.reset(static_cast<BasicWind*>(new testCq0Shields1986Wind(*this)));
    } else if (args().disk->wind == "Janiuk2015") {
        wind_.reset(static_cast<BasicWind*>(new Janiuk2015Wind(*this)));
    } else if (args().disk->wind == "Shields1986") {
        wind_.reset(static_cast<BasicWind*>(new Shields1986Wind(*this)));
    } else if (args().disk->wind == "Woods1996AGN") {
        wind_.reset(static_cast<BasicWind*>(new Woods1996AGNWind(*this)));
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
	set_Mdot_outer_boundary(obtain_Mdot_outer_boundary());
	verify_disc_mass(tau);
	invalidate_optional_structure();
	
	current_.i_t ++;
	current_.t += tau;
	wind_->update(*this);
}


void FreddiState::verify_disc_mass (double tau) {
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES-5) {std::cout << "+++ Mdisk("<< current_.t<<") = " << Mdisk() << " Rlast+1=" << R()[last()+1] << " delta_M_shift= " <<2.0 * M_PI * R()[last()] * Sigma()[last()] *(R()[last()+1]-R()[last()]) << "\n" << std::endl;               }
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES-5) {std::cout << "+++2 Mdot_in=" << Mdot_in() << " Mdot_out=" << obtain_Mdot_outer_boundary()<<" Mdisk+Mdot*tau=" << Mdisk()+(Mdot_in() - obtain_Mdot_outer_boundary()) * tau << "\n" << std::endl;               }
	if ( (args().calc->verb_level > VERB_LEVEL_MESSAGES-5) && (R()[last()+1]>0)) {std::cout << "+++  Mdisk+Mdot*tau+Sigma*2*Pi*R*dR =" << Mdisk()+(Mdot_in() + Mdot_in_prev()*2.3)  * tau  + 2.0 * M_PI * R()[last()] * Sigma()[last()] *(1.-3./4. * (R()[last()+1]-R()[last()])/(R()[last()+1]+R()[last()])  ) * (R()[last()+1]-R()[last()]) << "\n+++" << (1.-3./4. * (R()[last()+1]-R()[last()])/(R()[last()+1]+R()[last()])  ) * (R()[last()+1]-R()[last()]) << std::endl;               }
}



double FreddiState::obtain_Mdot_outer_boundary() {
    
    if (args().disk->DIM_front_approach == "outflow") {
	//if (args().disk->boundcond == "Hameury") {
	if (args().disk->scatter_by_corona == "yes") {
	    // assuming there is scatter above the disc
	    if (last() == Nx()-1) {
		return Mdot_out();
	    } else  {
		return -1.0*args().disk->DIM_front_Mdot_factor*Mdot_in();
		//-1.5 make Rfront/RMdot0 = 2
		// -2.5; -2.3 makes very close to Hameury with =Tcrit = (8840. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ; and Tcrit = 10300.;
	    }
	}
	//if (args().disk->boundcond == "Hameury_no_scatter") {
	if (args().disk->scatter_by_corona == "no") {
	    // if there is no scattering, the disc zone beyond maximum of Fvis (where dotM=0) is in the shadow
	    // This means its temperature
	    // if irradiation temperature is greater than critical, disc cannot be cold
	    // Tirr is checked at Thot or Tfront, depending on option boundcond (scatter/no scatter),see Tirr_critical. 
	    // If there is no scattering, the disc beyond r, where dotM=0, is in the shadow from the direct central radiation
	    
	    if ((last() == Nx()-1) or  (  Tirr().at(last()) >= Tirr_critical (R().at(last()), last()) ) ) {
		return Mdot_out();
	    } else  {
	   //     std::cout << "CF!" << std::endl;
		return -1.0*args().disk->DIM_front_Mdot_factor*Mdot_in();
		//-1.5 make Rfront/RMdot0 = 2
		// -2.5; -2.3 makes very close to Hameury with =Tcrit = (8840. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ; and Tcrit = 10300.;
	    }
	}
	
    } 
    // DIM_front_approach == "maxFvis" :
    return Mdot_out();
    
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
		const vecd& Shad = Shadow();
		const double Lbol = Lbol_disk();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = (1.0 - Shad[i]) * K[i] * Lbol * angular_dist_disk(H[i] / R()[i]) / (4. * M_PI * m::pow<2>(R()[i]));
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


const vecd& FreddiState::Shadow() {
        vecd x(Nx());
        double max_H2R = 0.0;
        for (size_t i = first(); i < Nx(); i++) {
        	const double Hrelative= Height()[i]/R()[i];
        	max_H2R = std::max(Hrelative, max_H2R);
		if (Hrelative >= max_H2R) {
       			x[i] = 0.0;
		} else {
		    if (args().disk->scatter_by_corona == "no") {
			x[i] = 1.0;
		    }
		}
		// x[i] = 0.0; // @XRPCALCApril24
        }
        //std::cout << "max= " << max_H2R << "H2R last = " << Height()[last()]/R()[last()]  << "last = " << last() << std::endl;
        opt_str_.Shadow = std::move(x);
        return *opt_str_.Shadow;
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
		C_w(state.args().disk->windparams.at("C_w")),
		R_w(state.args().disk->windparams.at("R_w")) {}

void FreddiState::testCq0Shields1986Wind::update(const FreddiState& state) {
	BasicWind::update(state);
	const double h_wind_min = std::sqrt(R_w) * state.h().back();
	for (size_t i = state.first(); i <= state.last(); ++i) {
		if (state.h()[i] > h_wind_min) {
			C_[i] = -0.5/M_PI * C_w * state.Mdot_in() /
					(std::log(1 / R_w) * m::pow<2>(state.R()[i])) *
					(4 * M_PI * m::pow<3>(state.h()[i])) / m::pow<2>(state.GM());
		}
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
       //C_[i] =  B_1 * 3000.0 * v_r * R_g * fout * ( std::pow((state.R()[state.last()]/R_g), 0.2) - std::pow((state.R()[i]/R_g), 0.2) );
       //*( std::pow((state.R()[state.last()]/R_g), 0.2) - std::pow((state.R()[i]/R_g), 0.2)
            const double fout = 1 - 1/(1+A_0*m::pow<2>(lol));
            const double C_0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
            const double Q = 2 * (3/(8 * M_PI)) * ((m::pow<4>(state.GM()))/(m::pow<7>(state.h()[i])))  ;
            B_[i] =  - 0.75 * C_0 * (1/B_1) * ((4 * Q * state.R()[i])/(3* state.GM()))*fout ; }
    }
}


FreddiState::Shields1986Wind::Shields1986Wind(const FreddiState& state):
        BasicWind(state),
        Xi_max(state.args().disk->windparams.at("Xi_max")),
        T_ic(state.args().disk->windparams.at("T_ic")),
        Pow(state.args().disk->windparams.at("Pow")) {
    update(state);
}

void FreddiState::Shields1986Wind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    //  1983ApJ...271...70B page 4
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    //  1983ApJ...271...70B page 3
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_ic);
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double L_crit = (1.0 / 8.0) * std::sqrt(GSL_CONST_CGSM_MASS_ELECTRON / (disk->mu * GSL_CONST_CGSM_MASS_PROTON)) * std::sqrt((GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT ) / (GSL_CONST_CGSM_BOLTZMANN * T_ic)) * L_edd;
    const double P_iC = L / (4.0 * M_PI * m::pow<2>(R_iC) * Xi_max * GSL_CONST_CGSM_SPEED_OF_LIGHT);
    const double el = L/L_crit;
    
    	std::cerr << R_iC << "\t" << disk->mu << "\t" << P_iC << "\t" << state.eta() << "\t" << state.Mdot_in() << "\t" <<  L << "\t" << L_edd << "\t" << L_crit << "\t" <<  GSL_CONST_CGSM_MASS_ELECTRON << std::endl; 
	//std::cerr << "\t" << L << "\t" << R_iC  << "\t" << L_edd << "\t" << C_ch << "\t" << P_0 << "\t" << m_ch0 << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        //  1986ApJ...306...90S page 2
        if (state.R()[i] > 0.1*R_iC) {
            const double xi = state.R()[i] / R_iC;
            const double T_ch = T_ic * std::pow(el, 2.0 / 3.0) * std::pow(xi, -2.0 / 3.0);
            const double P_0 = L / (4.0 * M_PI * m::pow<2>(state.R()[i]) * Xi_max * GSL_CONST_CGSM_SPEED_OF_LIGHT);
            const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ch) / (disk->mu * GSL_CONST_CGSM_MASS_PROTON));
            const double m_ch0 = P_0/C_ch;
            const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
            //  1986ApJ...306...90S appendix B page 16
            const double y = std::sqrt(1. + 1./(4.*m::pow<2>(xi)) + ((m::pow<2>(xi))/(1. + m::pow<2>(xi)))*((1.2*xi/(xi + el) + 2.2/(1. + m::pow<2>(el)*xi))*(1.2*xi/(xi + el) + 2.2/(1. + m::pow<2>(el)*xi))));
            const double Mach_cc = std::cbrt(((1. + (el + 1.)/xi)/(1. + 1./((1. + m::pow<2>(xi))*m::pow<4>(el)))));
	    const double p_po = 0.5 * std::exp(-(((1. - 1./y)*(1. - 1./y))/(2.*xi)));
	    C_[i] = - 2.0 * Pow * C0 * m_ch0 * Mach_cc * p_po * m::pow<2>(y);		}
    }
}


FreddiState::Woods1996AGNWind::Woods1996AGNWind(const FreddiState& state):
        BasicWind(state),
        C_0(state.args().disk->windparams.at("C_0")),
        T_ic(state.args().disk->windparams.at("T_ic")) {
    update(state);
}

void FreddiState::Woods1996AGNWind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_ic);
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
        T_ic(state.args().disk->windparams.at("T_ic")),
        Pow(state.args().disk->windparams.at("Pow")) {
    update(state);
}

void FreddiState::Woods1996ShieldsApproxWind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double R_iC = (state.GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON)/(GSL_CONST_CGSM_BOLTZMANN * T_ic);
    //const double VeL = std::sqrt(state.GM()/R_iC) ;
    //const double C_iC = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ic)/( GSL_CONST_CGSM_MASS_PROTON));
    //const double m_ch0 = disk->Mdot0 / (M_PI * state.R().back()*state.R().back());
    const double L_edd = (4.0 * M_PI * state.GM()* 2.0 * disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
    const double L_crit = (1.0 / 8.0) * std::sqrt(GSL_CONST_CGSM_MASS_ELECTRON / (disk->mu * GSL_CONST_CGSM_MASS_PROTON)) * std::sqrt((GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT* GSL_CONST_CGSM_SPEED_OF_LIGHT ) / (GSL_CONST_CGSM_BOLTZMANN * T_ic)) * L_edd;
    const double el = L/L_crit;
    
   //std::cerr << "\t" << Lcrit << "\t" << L_crit   << "\t" << std::endl;

    //std::cerr << "\t" << L << " \t" << L_edd << "\t" << L/L_edd  << "\t" << L_crit << "\t" << R_iC << "\t" << state.eta() << "\t" << std::endl;

    for (size_t i = state.first(); i <= state.last(); ++i) {
        if (state.R()[i] > 0.1*R_iC) {
            //  1986ApJ...306...90S page 2
            const double xi = state.R()[i] / R_iC;
            const double xi1 = R_iC / state.R()[i];
            const double T_ch = T_ic * std::pow(el, 2.0 / 3.0) * std::pow(xi, -2.0 / 3.0);
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
            C_[i] = - 2.0 * Pow * C0 * Fr * Fc * Expo;
        }

    }
}

FreddiState::PeriodPaperWind::PeriodPaperWind(const FreddiState& state):
	BasicWind(state),
	C_w(state.args().disk->windparams.at("C_w")) {
    update(state);
}

void FreddiState::PeriodPaperWind::update(const FreddiState& state) {
    BasicWind::update(state);
    const auto disk = state.args().disk;

    for (size_t i = state.first(); i <= state.last(); ++i) {
	//const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
	const double Mdot = state.Mdot_in() * (state.h()[i] - state.h()[state.first()]) /(m::pow<2>(state.h()[state.last()] - state.h()[state.first()])) ;
	C_[i] = - 2.0 * C_w * Mdot;
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
	// maximum density on the cold branch
	// Lasota et al., A&A 486, 523–528 (2008), Eq A.1, DOI: 10.1051/0004-6361:200809658
	return 74.6 * std::pow(args().basic->alphacold / 0.1, -0.83) * std::pow(r / 1e10, 1.18)
		* std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.40);
}

double FreddiState::Sigma_plus(double r) const { 
        // minimum density on the hot branch
	// Lasota et al., A&A 486, 523–528 (2008), Eq A.1, DOI: 10.1051/0004-6361:200809658
	return 39.9 * std::pow(args().basic->alpha / 0.1, -0.80) * std::pow(r / 1e10, 1.11)
		* std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.37);
}

//@@@@@ 28/09/2023
double FreddiState::Teff_plus(double r) const {
    // minimum temperature on the hot branch; Tavleev+23:
    double  A = 7341.;
    double beta = 0.0290;
    double gamma = -0.00484;
    double delta = -0.08426;
    return A * std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, beta) * std::pow(args().basic->alpha , gamma) * std::pow(r / 1e10, delta);
    
    // Lasota et al., A&A 486, 523–528 (2008), Eq A.1, DOI: 10.1051/0004-6361:200809658
    return 6890 *  std::pow(r / 1e10, -0.09) 
    * std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, 0.03);
}

// double FreddiState::Teff_minus(double r) const {
//     // maximum temperature on the cold branch Tavleev+23:
//     double  A = 6152.;
//     double beta = 0.0315;
//     double gamma = 0.00165;
//     double delta = -0.08977;
//     return A * std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, beta) * std::pow(args().basic->alpha , gamma) * std::pow(r / 1e10, delta);
//     
// }

double FreddiState::v_visc(double r, double z2r_at_r) {
	//oprel().Height(R()[i], F()[i])
	return  args().basic->alpha/oprel().Pi3 *  m::pow<2>(z2r_at_r)  * std::sqrt(GM()/r); 
	//Height2R()[]
}

double FreddiState::R_vis_struct (double r, double z2r_at_r)  {
        // estimate the radius, to which viscous process can reach in time tau
        return  r - v_visc(r, z2r_at_r) * args().calc->tau;       
}


double FreddiState::v_cooling_front(double r, double sigma_at_r) {
        // The cooling-front velocity depends on the ratio between the current Sigma and critical Sigmas
        // Ludwig et al., A&A 290, 473-486 (1994), section 3
        // https://ui.adsabs.harvard.edu/abs/1994A%26A...290..473L
        // units: cm/s
        const double Sigma_plus_ = Sigma_plus(r);
        
        //double sigma =  std::log( Sigma()[last()] / Sigma_plus_ ) /  std::log( Sigma_minus(r)/Sigma_plus_ ) ;
        //sigma = 0; 
	//Sigma_plus is the minimum density on the hot branch
        
        double sigma =   std::log( sigma_at_r / Sigma_plus_ ) /  std::log( Sigma_minus(r)/Sigma_plus_ ) ;
	if (sigma<0.) {sigma=0.0;}

        double v = 1e5 * (1.439-5.305*sigma+10.440*m::pow<2>(sigma)-10.55*m::pow<3>(sigma)+4.142*m::pow<4>(sigma))
               * std::pow(args().basic->alpha / 0.2, 0.85-0.69*sigma) 
               * std::pow(args().basic->alphacold / 0.05, 0.05+0.69*sigma)
               * std::pow(r / 1e10, 0.035)
               * std::pow(args().basic->Mx / GSL_CONST_CGSM_SOLAR_MASS, -0.012);
        if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_R v = " << v << " r=" << r << " Sigma="<< sigma_at_r << " sigma=" << sigma<< "\n" << std::endl;               }
	
        return v;
}

double FreddiState::R_cooling_front (double r, double sigma_at_r)  {
        // estimate the radius, where cooling front could reach, starting from r, if it moved with corresponding velocity
        // previous location of Rhot moves with the cooling-front velocity:
        return  r - v_cooling_front(r, sigma_at_r) * args().calc->tau;       
	
	
        return  R()[last()] - v_cooling_front(r,sigma_at_r) * args().calc->tau;   
        //return  R()[last()] - v_cooling_front(R()[last()]) * args().calc->tau  ; 
        // this variant leads to more abrupt evolution, since the front velocity is larger
}



double FreddiState::Rfront_Rhot (double r, double z_r) const {

	if (args().disk->Rfront_Mdotzero_factor == 1.0) {return 1.0;}
	// call: Rfront_Rhot(R().at(ii-1)) or Rfront_Rhot( R().at(ii) )
	double Rfront_Rhot=1.0;
	double Rfront_Rhot_variable = std::min(args().disk->Rfront_Mdotzero_factor, args().basic->rout/r);
        
       /* if (Rfront_Rhot_variable < args().disk->Rfront_Mdotzero_factor) { 		       
	    if (current_.R_dotM0_before_shift == 0.) {
		set_R_dotM0_before_shift(r);
	}
       } */
	/*	
		// at start: (1) R_dotM0 = Rout or (2) Rhot_prev < Rout, (2C)  if 2 Rhot_prev > Rout
		// (1-) no scatter -> R_dotM0(t) moves with v_visc;
// 		      (a) Qirr Yes : Rfront = R_dotM0(t)
// 		      (b) Qirr No : * Rfront increases from R_fotM0(t) until factor = 2
		// (1+) scatter ->  
// 		 	(a) Qirr Yes : * Rfront increases from R_fotM0(t) until factor = 2
// 		        (b) Qirr No : * Rfront increases from R_fotM0(t) until factor = 2
// 		(2C) At the beginning disc should collapse:
// 		(2C-) no scatter :
// 		    (a) Qirr Yes : Rhot = R_dotM0
// 		    (b) Qirr No : Rhot <= 2  R_dotM0 at t=0 and optionally still moves
// 		(2C+) scatter:
// 		    Qirr Yes : * Rhot <= 2 R_dotM0  at t=0 and optionally still moves
// 		    Qirr No : * Rhot <= 2 R_dotM0  at t=0 and optionally still moves
// 		(2-) no scatter 
// 		    Qirr Yes :
// 		    Qirr No :
// 		(2+) scatter
// 		    Qirr Yes :
// 		    Qirr No :
		// Fmaximum cannot move faster than v_visc
		// notice initial position of Fmaximum
		return Rfront_Rhot_variable;
		
	    }
	    if (current_.R_dotM0_before_shift > 0.) && (r > current_.R_dotM0_before_shift/args().disk->Rfront_Mdotzero_factor) {
		// Fmaximum moves until it shifts so that Rfront_Rhot_variable = 2:
		double R_dotM0 = R_vis_struct(current_.R_dotM0_before_shift,z_r)
		
		Rfront_Rhot_variable = std::max (args().disk->Rfront_Mdotzero_factor, current_.R_dotM0_before_shift/R_dotM0);
	    }
	}
	*/
	//if (args().disk->boundcond == "scatter_by_corona") {
	if (args().disk->scatter_by_corona == "yes") {    
	    Rfront_Rhot = Rfront_Rhot_variable;
	} 
	//if (args().disk->boundcond == "no_scatter_by_corona") {
	if (args().disk->scatter_by_corona == "no") {    
	    if (current_.maxR_Qirr_no_role == 0.) {
		// no scattering and irradiation's role is significant
		Rfront_Rhot = 1.; // the disk beyond dot M = 0 is shadowed 
	    } else {
		// no scattering and irradiation's role is NOT significant
		// then conditions of the hot zone are checked  at radius > r (dotM=0)
		Rfront_Rhot_variable = std::min( args().disk->Rfront_Mdotzero_factor, current_.maxR_Qirr_no_role / r);
		Rfront_Rhot = Rfront_Rhot_variable;
		//Rfront_Rhot=1.; //??????????????????
		
		// здесь надо написать так, чтобы 
	/*	if (Rfront_Rhot_variable < args().disk->Rfront_Mdotzero_factor) { 		       
		    if (current_.R_dotM0_before_shift == 0.) {
			set_R_dotM0_before_shift(r);
			// Fmaximum cannot move faster than v_visc
			// notice initial position of Fmaximum
			return Rfront_Rhot_variable;
			
		    }
		    if ((current_.R_dotM0_before_shift > 0.) && (r > current_.R_dotM0_before_shift/args().disk->Rfront_Mdotzero_factor)) {
			// Fmaximum moves until it shifts so that Rfront_Rhot_variable = 2:
			double R_dotM0 = R_vis_struct(current_.R_dotM0_before_shift,z_r)
			return std::max (args().disk->Rfront_Mdotzero_factor, current_.R_dotM0_before_shift/R_dotM0);
		    }
		}
		*/
	    }
	}
	return Rfront_Rhot;
}


double FreddiState::Tirr_critical (double r, int ii)  {
    // determine critical level Tirr which keeps the ring hot.
    // if we analyse r=Rfront,  we should take into account that Rfront
    // is farther than   r ==  R().at(ii) 
    //   
    // it might be possible that function can take no arguments, if it calculates only at R()[last()] and use Tirr().at(last()) / Tph_vis().at(last())
    // newest approach: radius_popravka = 1
    
    if (args().disk->boundcond == "Teff") {
	// never gets here
	throw std::invalid_argument("Why do you ask critical irradiation temperature if you set boundcond == Teff? Choose boundcond among Tirr/no_scatter_by_corona/scatter_by_corona\n");
    }

    // henceforth, only two possibilities remain: boundcond = "Tirr" or "no_scatter_by_corona"
    
    double radius_popravka = Rfront_Rhot( r,  oprel().Height(R()[ii], F()[ii])/r);
    
    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout <<  " c_E  R_last=" << R()[last()] << " R_ii="<< R()[ii] << " R_dotM0_before_shift=" << current_.R_dotM0_before_shift << " pop ="<<  radius_popravka  <<"\n" << std::endl;}
    
    
    // Qirr/Qvis increases linearly with radius, if Qvis ~ r^(-3/4) which is not exactly true beyond radius where dotM=0, because Fvis is not proportional to h there
    //thus, generally, Qirr_Qvis is the lower estimate 
    double Qirr_Qvis = pow(Tirr().at(ii) / Tph_vis().at(ii),4.)*radius_popravka ; 
    
    // determine radius, inside which irradiation is not significant:
    if ( (Qirr_Qvis < pow(args().disk->Tirr2Tvishot,4.) ) && (current_.maxR_Qirr_no_role == 0.) ) { 
	set_maxR_Qirr_no_role (r); 
    }  
    
    double Tcrit; // value to return
    
    if (args().disk->check_Temp_approach == "const") { 
	 // since we calculate critical level for irradiation temperature, we use the fact that Tirr ~ R(-1/2) (This is abs true only if Cirr = const)
	Tcrit = args().disk->Thot * pow(radius_popravka,0.5);
	
// 	if (args().disk->boundcond == "no_scatter_by_corona") {
	if (args().disk->scatter_by_corona == "no_") {    // @XRPCALCApril24
	    // do not take into account Tirr if Qirr/Qvis< critical_value
	    // and override previous assignment; instant return
	    if (Qirr_Qvis < pow(args().disk->Tirr2Tvishot,4.)) {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_Y Qirr/Qvis<crit \n" << std::endl;}
		// like the critical temperature is always too big 
		return 1e10;
	    } else
		// ( Tirr )
	    {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_W Qirr/Qvis=\n" <<Qirr_Qvis <<"\n"<<  std::endl;} 
	    }
	}
	
    } else if ((args().disk->check_Temp_approach == "Tavleev")  || (args().disk->check_Temp_approach == "Hameury") ) { 
        
        //if (args().disk->boundcond == "no_scatter_by_corona") {
	if (args().disk->scatter_by_corona == "no_") { //   @XRPCALCApril24
	    // if there is no scattering, the Rfront is in the shadow
	    // do not take into account Tirr if Qirr/Qvis < critical_value:
	    if (Qirr_Qvis < pow(args().disk->Tirr2Tvishot,4.)) {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_V Qirr/Qvis=\n" <<Qirr_Qvis <<"\n"<<  std::endl;} 
		// cooling wave will be going with irradiation having no effect; 
		// return a large value of Tcrit
		// like the critical temperature is always too big 
		return 1e10;
	    }
	}

	if (Qirr_Qvis > 1.) { 
	        // since we calculate critical level for irradiation temperature, we use the fact that Tirr ~ R(-1/2) (This is completely true only if Cirr = const)
		Tcrit = (9040. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ;
		if (args().disk->check_Temp_approach == "Hameury") {
		    Tcrit = (8640. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ; 
		}
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_F Qirr/Qvis>1 with Tcrit(Rhot)=" << Tcrit << " ii="<<ii<<"\n" << std::endl;}

	} else {
		Tcrit = (9040. - 2216.) * pow(radius_popravka,0.5) ;
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_G Qirr/Qvis<1 with Tcrit(Rhot)=" << Tcrit << " ii="<<ii<<"\n" << std::endl;}

	}
	if (args().disk->check_Temp_approach == "Hameury") {
	    if (ii < Nx()-1) {
		// 
		// if Rhot is less than Rtid, then critical temperature is fixed:
		// this behaviour is established for solution to agree with Hameury'solution
		Tcrit =  10300.; // Hameury email
	    }
	}
    } else {
	throw std::invalid_argument("Wrong check_Temp_approac: Choose const or Tavleev or Hameury\n");
    }
    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_H pop ="<<  radius_popravka << " Tirr(Rhot) = " << Tirr().at(ii) << " and Tcrit=" << Tcrit<< "\n" << std::endl;}
    
    return Tcrit;
    
}

//int FreddiState::ring_state_vertical(const int ii) {
int FreddiState::check_ring_is_cold(const int ii) {    
    // returns 1 for hot, 0 for cold
    
    // R().at(ii), the radius where Mdot = 0, or ?
    // multiplied by radius_popravka, gives the hot zone radius
    // Sigma_minus is the maximum density on the cold branch  
    // Sigma_plus is the minimum density on the hot branch     = Sigma_min(Menou+1999)

    
    if (args().disk->boundcond == "Teff") {
	if (Tph().at(ii) >= args().disk->Thot) {
	    return 0; // HOT
	} else {
	    return 1; // COLD
	}
    }
    
    if ((args().disk->boundcond == "Tirr") || (args().disk->boundcond == "Tirr_Ham??eury")) {
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "boundcond=Tirr="<<Tirr_critical (R().at(ii), ii) << std::endl;}
	if (Tirr().at(ii) >= Tirr_critical (R().at(ii), ii)) {
	    set_R_dotM0_before_shift( R().at(ii) );
	    return 0; // HOT
	} else {
	    return 1; // COLD

	}
    }
   
   // only no_scatter_by_corona is possible ; everything else should be removed here:
    if (!(args().disk->boundcond == "DIM")) { 
	throw std::invalid_argument("Wrong boundcond at check_ring_is_cold() in freddi_state");
    }
    
    
    if (  Tirr().at(ii) >= Tirr_critical (R().at(ii), ii) ) {
	// if irradiation temperature is greater than critical, disc cannot be cold
	// Tirr is checked at Thot or Tfront, depending on option boundcond (scatter/no scatter),see Tirr_critical. 
	// If there is no scattering, the disc beyond r, where dotM=0, is in the shadow from the direct central radiation
	// In practice, temperature of the disc at R, where dotM=0, is checked against the critical temperature multiplied by some factor
	if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_K Tirr high - HOT \n" << std::endl;}
	// memorise radius before shift of the hot boundary
	set_R_dotM0_before_shift( R().at(ii) );
	return 0; // HOT
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
	    return 0; // HOT
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
		    return 1; // COLD
		} else {
		    //throw std::invalid_argument("check_ring_is_cold: logic mistake 1");
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
	    } 
// 	    else {
		
		//  R_cooling_front = r - v_cooling_front(r, sigma_at_r) * args().calc->tau;   
//   		if (radius_popravka <= args().disk->Rfront_Mdotzero_factor) {
// 		    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_PP pop ="<< radius_popravka <<  " radius_popravka < Rfront_Mdotzero_factor->COLD \n" << std::endl;}
//   		    return 1;
//   		}

		// check cooling front position; cooling front is moving from last Rout
		
		if  (current_.R_dotM0_before_shift == 0 ) {
		     // COLD ; just collapse to initial position:
		    return 1;
		}
		if ((radius_popravka * R().at(ii) > R_cooling_front ( Rfront_Rhot(R()[last()], oprel().Height(R()[last()], F()[last()])/R()[last()] ) * R()[last()] , Sigma().at(last())*sigma_factor) )  && ( last() < args().calc->Nx-1)  ){
		    //radius is beyond front 
		    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "++c_P pop ="<< radius_popravka <<  " beyond front - COLD \n" << std::endl;}
		    // COLD
		    return 1;
		} else {
		    if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "++c_Q pop ="<<  radius_popravka <<  " inwards front - HOT \n" << std::endl;}
		    // in fact front starts from a greater distance?
		    set_R_dotM0_before_shift( R().at(ii) );
		    // HOT
		    return 0;
		}
//	   }
	}
    }
    throw std::invalid_argument("check_ring_is_cold: logic mistake 2");
    return 0;
}

#undef VERB_LEVEL_MESSAGES


