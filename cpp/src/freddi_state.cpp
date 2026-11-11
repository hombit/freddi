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
		v_esc(initialize_v_esc(h, GM)),
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

vecd FreddiState::DiskStructure::initialize_v_esc(const vecd& h, double GM) {
    vecd v_esc(h.size());
    for (size_t i = 0; i < h.size(); i++) {
        // v_esc = sqrt( 2 * G * M / R )
		v_esc[i] = std::sqrt(2.0) * GM / h[i];
    }
    return v_esc;
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



double FreddiState::obtain_Mdot_outer_boundary() const {
    
    
    if (last() == Nx()-1) {
	// hot radius is at the outer geometrical radius
	// Rhot = Rtid
	return Mdot_out();
    }
    
    // if the hot radius is less than geometrical and there is a cold ring, 
    // Rhot<Rtid
    // then it is important whether there is a scattering corona or not
    //                      and if flag DIM_front_approach == "outflow"
  
    if (args().disk->DIM_front_approach == "outflow") {
	// there is an outward mass flow at the front; this negative rate  = factor * inner accretion rate during self-similar stage
	
		if (args().disk->scatter_by_corona == "yes") {
			// there is scattering above the disc
			// there is actually an option to add Mdot_out() here .....
			return -1.0*args().disk->DIM_front_Mdot_factor*Mdot_in();
			//-1.5 make Rfront/RMdot0 = 2
			// -2.5; -2.3 makes very close to Hameury with =Tcrit = (8840. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ; and Tcrit = 10300.;
		}
		
		if (args().disk->scatter_by_corona == "no") {
			// if there is no scattering, the disc zone beyond maximum of Fvis (where dotM=0) is in the shadow from the direct central radiation
			// 
			// if irradiation temperature at Rhot is greater than critical, disc cannot be cold; but the ring outside Fvis=max is cold
			// that means  that  Mdot = 0  at Rhot;
			// Tirr is checked at Rhot,see Tirr_critical. 
			
			if (  Tirr().at(last()) >= Tirr_critical (R().at(last()), last())  ) {
			return (0.);
			} else  {
			// if irradiation is unable to keep the disc hot, there is a smooth transition zone up to R hot and there is an outflow at Rhot
			return -1.0*args().disk->DIM_front_Mdot_factor*Mdot_in();
			//-1.5 make Rfront/RMdot0 = 2
			// -2.5; -2.3 makes very close to Hameury with =Tcrit = (8840. - 2216.* 1./Qirr_Qvis ) * pow(radius_popravka,0.5) ; and Tcrit = 10300.;
			}
		}
	
    } 
    
    // if DIM_front_approach is not DIM, then just return transfer rate
    // DIM_front_approach == "maxFvis" : Mdot = 0 at the front; it is actually is set as incoming to the disc Mdot_ot;
    return Mdot_out();
    
}

double FreddiState::Mdot_in() const {
	//std::cout << "****\n" << std::endl;
	return (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
}


double FreddiState::Mdot_out_from_F() const {
	//std::cout << "*calc out Mdot ***\n" << std::endl;
	return (F()[last()] - F()[last()-1] ) / (h()[last()] - h()[last()-1]);
}

double FreddiState::Lbol_disk() const {
	return eta() * Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
}


double FreddiState::L_edd_disk() const {
	return 4.0 * M_PI * GM() * 2.0 * args().disk->mu * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT /
			GSL_CONST_CGSM_THOMSON_CROSS_SECTION;
}


double FreddiState::T_ic_current(double T_ic_const) const {
	// windT_ic_approach: how the inverse Compton temperature of the thermal winds is determined
	if (args().disk->windT_ic_approach == "const") {
		return T_ic_const;
	}
	if (args().disk->windT_ic_approach == "Done2018") {
		//  Done, Tomaru & Takahashi 2018 (MNRAS 473, 838), Eqs. (8)-(9), smoothly blended
		//  around l = L/Ledd = 0.02 following Dubus et al. (2019), Eq. after their Eq. (9):
		//  T_IC = T_HS/(1+x) + T_SS*x/(1+x), x = (l/0.02)^6
		const double l = Lbol_disk() / L_edd_disk();
		const double T_hs = (4.2 - 4.6 * std::log10(l / 0.02)) * 1e7;
		const double T_ss = 0.36 * std::pow(l / 0.02, 0.25) * 1e7;
		const double x = std::pow(l / 0.02, 6.0);
		return T_hs / (1.0 + x) + T_ss * (x / (1.0 + x));
	}
	throw std::invalid_argument("Wrong windT_ic_approach");
}


double FreddiState::R_iC_from_T_ic(double T_ic) const {
	double R_iC = GM() * args().disk->mu * GSL_CONST_CGSM_MASS_PROTON / (GSL_CONST_CGSM_BOLTZMANN * T_ic);
	if (args().disk->windT_ic_approach == "Done2018") {
		// Dubus et al. (2019) Eq. (4): radiation-pressure correction to R_IC, gated on
		// windT_ic_approach=Done2018 since that's the T_ic(l) relation this correction was derived
		// alongside; applied here (not just in the Cirr_Dubus2019()/Mdot_wind_Dubus2019()
		// diagnostics) so the actually simulated wind uses the same R_IC Dubus et al. (2019) do.
		// For L/Ledd >= 1/sqrt(2) the factor is <= 0 (radiation pressure unbinds the gas at all
		// radii, i.e. the wind is launched everywhere): clamp to a small positive floor
		// (1e-6 of the uncorrected value) so that xi = R/R_IC stays finite and positive.
		const double correction = 1.0 - std::sqrt(2.0) * Lbol_disk() / L_edd_disk();
		R_iC *= std::max(correction, 1e-6);
	}
	return R_iC;
}


double FreddiState::R_wind_inner_launch_radius(double R_iC) const {
	// Inner radius beyond which the thermal wind is launched, as a multiple of the passed-in
	// Compton radius R_iC. Begelman et al. (1983)/Shields et al. (1986)/Woods et al. (1996), and
	// Freddi's long-standing default, take this to be 0.1 R_IC; Dubus et al. (2019) instead assume
	// the wind is effective for R >= 0.2 R_IC. Controlled by --windR_launch_factor (windparams
	// "R_launch_factor"), default 0.1 -- set to 0.2 to compare against Dubus et al. (2019). Only
	// meaningful for --windtype=Woods1996 (reads its windparams).
	return args().disk->windparams.at("R_launch_factor") * R_iC;
}


double FreddiState::phase_opt() const {
	return 2.0 * M_PI * (t() - args().flux->ephemeris_t0) / args().basic->period;
}


double FreddiState::Lx() const {
	if (!opt_str_.Lx) {
		opt_str_.Lx = Luminosity(Tph_X(), args().flux->emin, args().flux->emax) / m::pow<4>(args().flux->colourfactor);
	}
	return *opt_str_.Lx;
}


double FreddiState::Fx() const {
    return Lx() * angular_dist_disk(cosi()) / (FOUR_M_PI * m::pow<2>(distance()));
}


const vecd& FreddiState::W() const {
	if (!opt_str_.W) {
		auto x = wunc()(h(), F(), first(), last());
		x.resize(Nx(), 0.0);
		opt_str_.W = std::move(x);
	}
	return *opt_str_.W;
}


const vecd& FreddiState::Sigma() const {
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


const vecd& FreddiState::Tph() const {
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


const vecd& FreddiState::Tirr() const {
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


const vecd& FreddiState::Qx() const {

	if (!opt_str_.Qx) {
		vecd x(Nx());
		const vecd& K = Kirr(); // Shadow is included in Kirr, see Kirr()
		const vecd& H = Height();
		const double Lbol = Lbol_disk();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = K[i] * Lbol * angular_dist_disk(H[i] / R()[i]) / (4. * M_PI * m::pow<2>(R()[i]));
		}
		opt_str_.Qx = std::move(x);
	}
	return *opt_str_.Qx;
}


double FreddiState::Cirr_scatter_dependent(size_t i) const {
	// uniform formula for the whole disc
	// Cirr = etaX * wind_Column_density *kappa_Thomson * angular_distribution / 2, see Qx()
	// here we calculate Kirr = etaX * wind_Column_density *kappa_Thomson / 2
	// scattering_opacity usually is kappa_Thomson
	//  $X=0.735$ $Y=0.248$ and $Z=0.017$, with $Z/X$ = 0.023 (GREVESSE and A.J. SAUVAL 1988)
	return args().irr->etaX * args().irr->scattering_opacity * Column_density_wind()[i] / 2;
}


double FreddiState::Cirr_Dubus2019() const {
	// Undocumented (comparison-only): Dubus (2019)-style wind-scattering irradiation
	// parameter, a single disc-wide constant instead of scatter_dependent's per-ring
	// Column_density_wind() estimate. Their Eq. (10):
	//   C = int_0^1 int_{Rin}^{Rout} sigma_T n_w mu dmu dr ~= kappa * Mdot_wind / (8 pi Rin v_w)
	// where kappa = sigma_T/m_I is the scattering opacity (--scattering_opacity); Rin = 0.2
	// R_IC is the wind launching radius (their Eq. 4, "thermal wind is effective for radii
	// R >= 0.2 R_IC"); Mdot_wind is the total (two-sided) wind mass loss rate, i.e. their
	// Mdot_w(Rin) integrated out to Rout, matching Freddi's existing Mdot_wind(); and v_w is
	// their "mass-weighted wind outflow rate" (text below their Eq. 10) -- the local escape
	// velocity v_esc(r) = sqrt(2GM/r) averaged (weighted by dMdot_wind/dr, not dMdot_wind_dh()'s
	// native per-dh density -- converted via dh/dr = GM/(2h), h = sqrt(GM r)) over the actual
	// radial wind mass-loss profile, not v_esc evaluated at the single radius Rin. The two can differ a lot:
	// the simulated wind's radial extent is set by --windR_launch_factor (only used by
	// --windtype=Woods1996, via R_wind_inner_launch_radius() -- not by this Rin), which need not
	// match their fixed 0.2 R_IC. E.g. for wind_Dubus_fig3.ini's windR_launch_factor=1e-3, Rin=0.2
	// R_IC there even exceeds Rout, so nearly all of Mdot_wind() actually originates well inside
	// Rin, and a v_w evaluated only at Rin would not represent where the wind mass is.
	// R_IC here is their Eq. (4): GM*mu*mp/(kB*T_ic) with an added rough radiation-pressure
	// correction (1 - sqrt(2) L/Ledd), applied here unconditionally so this diagnostic always
	// reflects the Dubus2019 R_IC regardless of --windT_ic_approach. (When
	// --windT_ic_approach=Done2018, R_iC_from_T_ic() already applies this same correction to the
	// actually simulated wind's own R_IC -- computed independently here via wind_->T_ic() rather
	// than wind_->R_IC() to avoid double-applying it in that case.) Per the source paper this
	// correction is explicitly "rough" and can go unphysical (negative) as L approaches
	// Ledd/sqrt(2) -- not guarded against here.
	const double R_IC_base = GM() * args().disk->mu * GSL_CONST_CGSM_MASS_PROTON / (GSL_CONST_CGSM_BOLTZMANN * wind_->T_ic());
	const double R_IC_corrected = R_IC_base * (1.0 - std::sqrt(2.0) * Lbol_disk() / L_edd_disk());
	const double R_in = 0.2 * R_IC_corrected;
	const double Mdot_w = Mdot_wind();
	if (Mdot_w <= 0.0) {
		return 0.0;
	}
	const vecd& Vesc = v_esc();
	const vecd& H = h();
	const double GM_ = GM();
	// dMdot_wind_dh() is a density per unit h; convert to a density per unit r via
	// dh/dr = GM/(2h) (from h = sqrt(GM r)) so the average below is mass-weighted in r.
	auto dMdot_wind_dr = [this, &H, GM_](size_t i) { return dMdot_wind_dh(i) * GM_ / (2.0 * H[i]); };
	const double Mdot_w_dr = trapz(R(), dMdot_wind_dr, first(), Nx() - 1);
	const double v_w = trapz(R(), [&Vesc, &dMdot_wind_dr](size_t i) { return Vesc[i] * dMdot_wind_dr(i); }, first(), Nx() - 1) / Mdot_w_dr;
	return args().irr->scattering_opacity * Mdot_w / (8.0 * M_PI * R_in * v_w);
}


double FreddiState::Cirr_direct_analytic(size_t i) const {
	// CIRR = eta_X * (dz/dr - z/r) * angular_distribution
	if (i == first()) {
		return 0.0; // no irradiation at the inner boundary
	}
	return (1.0 - Shadow()[i]) * args().irr->etaX *
			((Height()[i] - Height()[i - 1]) / (R()[i] - R()[i - 1]) - Height()[i] / R()[i]);
}


const vecd& FreddiState::Kirr() const {

	if(!opt_str_.Kirr) {
		vecd x(Nx());
		const vecd& Shad = Shadow();
		if (args().irr->irradiation_type == "scatter_dependent") {
			for (size_t i = first(); i < Nx(); i++) {
				x[i] = Cirr_scatter_dependent(i);
			}
		} else if (args().irr->irradiation_type == "Dubus2019") {
			const double C = Cirr_Dubus2019();
			for (size_t i = first(); i < Nx(); i++) {
				x[i] = C;
			}
		} else if  (args().irr->irradiation_type == "direct_analytic") {
			for (size_t i = first(); i < Nx(); i++) {
				x[i] = Cirr_direct_analytic(i);
			}
		} else if  (args().irr->irradiation_type == "constant_Cirr") {
			const vecd& H = Height();
			// if scatter_by_corona = yes, then there is no shadow
			// if there is scattering, then the cold part of the disc is irradiated by scattered radiation, so we use the same formula as for hot part
			for (size_t i = first(); i <= last(); i++) {
				x[i] =  args().irr->Cirr * std::pow(H[i] / (R()[i] * 0.05), args().irr->irrindex);
				if (args().disk->scatter_by_corona == "no") {
					x[i] *= (1.0 - Shad[i]);
				}
			}
			for (size_t i = last() + 1; i < Nx(); i++) {
				x[i] = args().irr->Cirr_cold * std::pow(H[i] / (R()[i] * 0.05), args().irr->irrindex_cold);
				if (args().disk->scatter_by_corona == "no") {
					x[i] *= (1.0 - Shad[i]);
				}
			}

		}
		opt_str_.Kirr = std::move(x);
	}
	return *opt_str_.Kirr;
}


const vecd& FreddiState::Height() const {
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


const vecd& FreddiState::Shadow() const {
        vecd x(Nx());
        double max_H2R = 0.0;
        for (size_t i = first(); i < Nx(); i++) {
        	const double Hrelative= Height()[i]/R()[i];
        	max_H2R = std::max(Hrelative, max_H2R);
			if (Hrelative >= max_H2R) {
       			x[i] = 0.0; // no shadow
			} else {
				x[i] = 1.0; // shadowed by the disc itself
			}
			// x[i] = 0.0; // @XRPCALCApril24
        }
        opt_str_.Shadow = std::move(x);
        return *opt_str_.Shadow;
}




const vecd& FreddiState::Tph_vis() const {
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

const vecd& FreddiState::Tph_X() const {
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


double FreddiState::lazy_magnitude(boost::optional<double>& m, double lambda, double F0) const {
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


double FreddiState::flux_star(const double lambda, const double phase) const {
	return star_.luminosity({inclination(), phase}, lambda) / (FOUR_M_PI * m::pow<2>(distance()));
}

double FreddiState::flux_star(const Passband& passband, const double phase) const {
	return star_.luminosity({inclination(), phase}, passband) / (FOUR_M_PI * m::pow<2>(distance()));
}


double FreddiState::dMdot_wind_dh(size_t i) const {
	double dFdh;
	if (i == first()) {
		dFdh = (F()[i+1] - F()[i]) / (h()[i+1] - h()[i]);
	} else if (i == Nx() - 1) {
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
}


double FreddiState::Mdot_wind() const {
	// Integrated over the whole disc (hot + cold), not just the hot region: Compton-wind mass loss
	// is driven by X-ray irradiation and R/R_IC, not by whether the local annulus is currently on
	// the hot or cold branch of the thermal instability, so a cold, quiescent outer disc can still
	// carry a wind if R >= windR_launch_factor*R_IC (Shields1986Wind, Woods1996AGNWind and
	// Woods1996ShieldsApproxWind all now fill windC() out to the outer edge accordingly). See also
	// Mdot_wind_cold() for the cold-zone-only contribution.
	if (!opt_str_.Mdot_wind) {
		opt_str_.Mdot_wind = trapz(h(), [this](size_t i) { return dMdot_wind_dh(i); }, first(), Nx() - 1);
	}
	return *opt_str_.Mdot_wind;
}


double FreddiState::Mdot_wind_cold() const {
	// The cold-zone-only (R > Rhot) contribution to Mdot_wind() above, exposed separately as a
	// diagnostic of how much of the total wind mass loss originates beyond the hot disc.
	if (!opt_str_.Mdot_wind_cold) {
		// trapz() returns 0 when first >= last, i.e. when there is no cold zone at all
		opt_str_.Mdot_wind_cold = trapz(h(), [this](size_t i) { return dMdot_wind_dh(i); }, last() + 1, Nx() - 1);
	}
	return *opt_str_.Mdot_wind_cold;
}


double FreddiState::Mdot_wind_Dubus2019() const {
	// Undocumented (comparison-only): total wind mass-loss rate under the Woods1996ShieldsApproxWind
	// formula (--windtype=Woods1996), with R_iC replaced by the Dubus (2019) radiation-pressure-
	// corrected Compton radius (their Eq. 4, same expression as in Cirr_Dubus2019()) instead of
	// Freddi's own R_iC = GM*mu*mp/(kB*T_ic). Purely diagnostic: does not feed back into the
	// evolution. The correction is applied unconditionally here (independent of
	// --windT_ic_approach) via wind_->T_ic() rather than wind_->R_IC(), since the latter may
	// already carry this same correction when --windT_ic_approach=Done2018 (see
	// R_iC_from_T_ic()) -- using it here too would double-apply it. Only meaningful for
	// --windtype=Woods1996 (reads its windparams); reports 0 (with a one-time warning) for any
	// other --windtype instead of crashing.
	const auto disk = args().disk;
	if (!disk->windparams.count("Xi_max") || !disk->windparams.count("Pow") ||
			!disk->windparams.count("IrAngDis") || !disk->windparams.count("R_launch_factor")) {
		static bool warned = false;
		if (!warned) {
			std::cerr << "Warning: Mdot_wind_Dubus2019 diagnostic output column is only meaningful "
					"for --windtype=Woods1996; reporting 0 since --windtype=" << disk->wind << std::endl;
			warned = true;
		}
		return 0.0;
	}
	const double Xi_max = disk->windparams.at("Xi_max");
	const double Pow = disk->windparams.at("Pow");
	const bool IrAngDis = disk->windparams.at("IrAngDis") != 0.0;
	const double T_ic = wind_->T_ic();
	const double L = Lbol_disk();
	const double L_edd = L_edd_disk();
	const double R_iC_base = GM() * disk->mu * GSL_CONST_CGSM_MASS_PROTON / (GSL_CONST_CGSM_BOLTZMANN * T_ic);
	const double R_iC = R_iC_base * (1.0 - std::sqrt(2.0) * L / L_edd); // Dubus (2019) Eq. (4)
	const double L_crit = (1.0 / 8.0) * std::sqrt(GSL_CONST_CGSM_MASS_ELECTRON / (disk->mu * GSL_CONST_CGSM_MASS_PROTON)) * std::sqrt((GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT* GSL_CONST_CGSM_SPEED_OF_LIGHT ) / (GSL_CONST_CGSM_BOLTZMANN * T_ic)) * L_edd;
	double el = L / L_crit;

	vecd C(Nx(), 0.0);
	// Whole disc, matching Mdot_wind()/Woods1996ShieldsApproxWind::update() -- see there for why.
	for (size_t i = first(); i < Nx(); ++i) {
		if (R()[i] > R_wind_inner_launch_radius(R_iC)) {
			if (IrAngDis) {
				el *= angular_dist_disk(Height()[i] / R()[i]);
			}
			const double xi = R()[i] / R_iC;
			const double xi1 = R_iC / R()[i];
			const double T_ch = T_ic * std::pow(el, 2.0 / 3.0) * std::pow(xi, -2.0 / 3.0);
			const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ch) / (disk->mu * GSL_CONST_CGSM_MASS_PROTON));
			const double C0 = (4.0 * M_PI * m::pow<3>(h()[i])) / (m::pow<2>(GM()));
			const double Fr = L / (4.0 * M_PI * m::pow<2>(R()[i]) * Xi_max * C_ch * GSL_CONST_CGSM_SPEED_OF_LIGHT);
			const double Fc = ((std::pow((1 + m::pow<2>(((0.125 * el + 0.00382) * xi1))), (1.0 / 6.0))) /
							   std::pow((1 + m::pow<-2>((m::pow<4>(el)* (1.0 + 262.0 * m::pow<2>(xi))))),
										(1.0 / 6.0)));
			const double Expo = std::exp(-(((1.0 - (1 / std::sqrt(1.0 + 0.25 * m::pow<2>(xi1)))) *
											(1.0 - (1 / std::sqrt(1.0 + 0.25 * m::pow<2>(xi1))))) / (2.0 * xi)));
			C[i] = -2.0 * Pow * C0 * Fr * Fc * Expo;
		}
	}

	auto dMdot_dh = [&C](const size_t i) -> double {
		return -C[i];
	};
	return trapz(h(), dMdot_dh, first(), Nx() - 1);
}


const vecd& FreddiState::Column_density_wind() const {
    if (!opt_str_.Column_density_wind) {
        vecd x(Nx(), 0.0);
        // Pre-calculate the running Mdot vector once
        const vecd& mdot_run = Mdot_wind_running();

        // Whole disc (not just the hot region), so the cold outer disc's own Cirr
        // (scatter_dependent) reflects any wind mass loss originating out there too.
        for (size_t i = first(); i < Nx(); i++) {
            // Formula: Σ = Mdot / (4 * pi * R * v_wind)
            x[i] = mdot_run[i] / (4.0 * M_PI * R()[i] * v_wind()[i]);
        }
        opt_str_.Column_density_wind = std::move(x);
    }
    return *opt_str_.Column_density_wind;
}

const vecd& FreddiState::tau_perp_wind() const {
    if (!opt_str_.tau_perp_wind) {
        vecd x(Nx(), 0.0);
        // tau_corona_perp = Column_density_wind * scattering_opacity, see the
        // --irradiation_type=scatter_dependent help text (options.cpp) for the derivation.
        const vecd& Sigma_wind = Column_density_wind();
        for (size_t i = first(); i < Nx(); i++) {
            x[i] = args().irr->scattering_opacity * Sigma_wind[i];
        }
        opt_str_.tau_perp_wind = std::move(x);
    }
    return *opt_str_.tau_perp_wind;
}

// Relative disc semi-height z_disk/R = h_d/xi assumed for the corona's source-to-disc optical
// depth tau(xi) (Ostriker, McKee & Klein 1991, Eqs. 39-40; see tau_r_corona()) in the cold zone
// beyond Rhot (Height() itself, which is irradiation-related, is not changed)
static constexpr double z_disk_to_R_cold_corona = 0.01;

double FreddiState::R_IC_corona() const {
    // Genuine Compton radius R_IC = GM mu m_p / (k_B T_ic) used for the Ostriker et al. (1991)
    // corona (xi = R/R_IC), WITHOUT the Dubus et al. (2019) Eq. (4) radiation-pressure
    // correction (1 - sqrt(2) L/Ledd), which R_iC_from_T_ic() / wind_->R_IC() apply to the wind
    // under --windT_ic_approach=Done2018. 0 for wind types without an inverse Compton temperature.
    const double T_ic = wind_ ? wind_->T_ic() : 0.0;
    if (T_ic <= 0.0) {
        return 0.0;
    }
    return GM() * args().disk->mu * GSL_CONST_CGSM_MASS_PROTON / (GSL_CONST_CGSM_BOLTZMANN * T_ic);
}

double FreddiState::R_hc() const {
    // Ostriker, McKee & Klein (1991, ApJ 377, 593) Eq. (20): the minimum radius R_hc where a hot
    // Compton-heated corona can exist, i.e. where the direct central flux L/(4 pi R^2) exceeds
    // twice the locally generated viscous flux 2 Qvis, with Qvis = 3 Fvis omega_K / (8 pi R^2)
    // (omega_K = sqrt(GM/R^3)). Qvis is evaluated from F() over the whole grid (not from Tph_vis,
    // which is only filled up to the hot-zone edge), so R_hc does not depend on Rhot. Both L and
    // Qvis scale with Mdot, so R_hc is nearly independent of Mdot (it depends on Mx, eta).
    // Qvis drops faster than R^-2 in the outer disc, so the direct flux wins outside R_hc: the
    // outermost crossing is searched for (log-linear interpolation between grid nodes).
    // Returns R()[first()] if the direct flux dominates everywhere, and +inf if nowhere.
    const double L = Lbol_disk();
    const size_t out = Nx() - 1;
    auto g = [&](size_t i) {
        const double omega_K = std::sqrt(GM() / m::pow<3>(R()[i]));
        const double Qvis = 3.0 * F()[i] * omega_K / (8.0 * M_PI * m::pow<2>(R()[i]));
        return L / (4.0 * M_PI * m::pow<2>(R()[i])) - 2.0 * Qvis;
    };
    if (g(out) <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    for (size_t i = out; i > first(); --i) {
        if (g(i - 1) <= 0.0) {
            const double g_in = g(i - 1), g_out = g(i);
            const double w = -g_in / (g_out - g_in);
            return std::exp((1.0 - w) * std::log(R()[i - 1]) + w * std::log(R()[i]));
        }
    }
    return R()[first()];
}

const vecd& FreddiState::tau_r_wind() const {
    // Optical (Thomson/scattering) depth along z = const, i.e. radially in the disc plane, through
    // the thermal wind from the wind launch radius to R:
    //   tau_r_wind(R) = scattering_opacity * int_{R_launch}^{R} rho_wind dR',
    //   rho_wind = Mdot_wind_running / (4 pi R^2 v_wind)  (so that int rho dR ~ Column_density_wind).
    // R_launch = R_wind_inner_launch_radius(R_IC) (windR_launch_factor * R_IC); for wind types
    // without it, integrated from the inner disc edge. Trapezoidal rule, starting from the first
    // grid node with R >= R_launch. Zero if R_IC <= 0.
    // NOTE: only the thermal (Compton-heated) wind is included. Other outflows
    // (super-Eddington winds, magnetic winds) are ignored.
    if (!opt_str_.tau_r_wind) {
        vecd x(Nx(), 0.0);
        const double R_IC = wind_ ? wind_->R_IC() : 0.0;
        if (R_IC > 0.0) {
            const double R_launch = args().disk->windparams.count("R_launch_factor")
                    ? R_wind_inner_launch_radius(R_IC) : 0.0;
            const vecd& mdot_run = Mdot_wind_running();
            const double kappa = args().irr->scattering_opacity;
            auto rho = [&](size_t i) {
                return mdot_run[i] / (4.0 * M_PI * m::pow<2>(R()[i]) * v_wind()[i]);
            };
            bool started = false;
            for (size_t i = first(); i < Nx(); i++) {
                if (R()[i] < R_launch) {
                    continue;
                }
                if (started) {
                    x[i] = x[i-1] + kappa * 0.5 * (rho(i) + rho(i-1)) * (R()[i] - R()[i-1]);
                }
                started = true;
            }
        }
        opt_str_.tau_r_wind = std::move(x);
    }
    return *opt_str_.tau_r_wind;
}

double FreddiState::z_disk_corona(size_t i) const {
    // Disc semi-height H_d = h_d R_IC used in the corona's source-to-disc optical depth tau(xi)
    // (OMK91 Eqs. 39-40, see tau_r_corona()): Height() in the hot zone (up to Rhot), and
    // z_disk_to_R_cold_corona * R in the cold zone beyond Rhot.
    return i <= last() ? Height()[i] : z_disk_to_R_cold_corona * R()[i];
}

const vecd& FreddiState::tau_r_corona() const {
    // Optical depth tau(xi) through the dynamically stable Compton-heated corona along the ray
    // from the central source to the disc surface at (xi, h_d(xi)), Ostriker, McKee & Klein
    // (1991, ApJ 377, 593; OMK91) Eqs. (32), (35). The ray is nearly parallel to the disc plane,
    // hence the column name "tau_r". This is NOT OMK91's tau_parallel (their Eqs. 48-49, the
    // depth from xi outward to the disc edge). It is an integral, not a local value:
    //   tau(xi) ~ int_{xi_min}^{xi} R_IC sigma_T n_e,b(xi') dxi',
    //   R_IC sigma_T n_e,b = (L/Ledd) / (Xi_tilde xi'^2)                          (Eqs. 29a, 37)
    //   => tau(xi) = (L/Ledd)/Xi_tilde * [1/xi_min - 1/xi]                        (Eq. 39)
    // where xi_min is where the straight ray enters the corona, i.e. crosses one coronal scale
    // height h_c(xi') = xi'^{3/2} (isothermal at T_ic, h_c/R = xi^{1/2}) above the base:
    //   (h_x - xi_min^{3/2}) / xi_min = (h_x - h_d(xi)) / xi                     (Eq. 40)
    // Here h_x = 0 (source in the disc plane; may change for the super-Eddington regime), for
    // which Eq. (40) is solved exactly: xi_min^{1/2} = h_d/xi = z_disk/R, so xi_min = (z_disk/R)^2
    // and Eq. (39) gives
    //   tau(xi) = (L/Ledd)/Xi_tilde * max(0, (R/z_disk)^2 - R_IC/R)              (exact, h_x = 0)
    // It is zero where xi <= (z_disk/R)^2, i.e. where the corona scale height is below the disc
    // semi-height (h_c < h_d): OMK's thick-corona assumption h_c >> h_d fails there and the ray
    // stays below the coronal base layer (isochoric model: no coronal density there).
    // Previously OMK's approximate solution of Eq. (40), Eq. (44), was used, i.e. Eq. (45):
    //   tau(xi) ~ (L/Ledd)/Xi_tilde / (h_x^{2/3} + (h_d/xi)^2)                   (Eq. 45)
    // which for h_x = 0 drops the -R_IC/R term; it overestimates tau by 1/(1 - xi_min/xi), and
    // spuriously gave tau ~ 1-3 in the inner disc where xi < xi_min.
    // WARNING: earlier versions of this comment and of the output help texts called the h_x = 0
    // form of Eq. (45) "(G.2)" and cited "Appendix G" of OMK91. OMK91 has only Appendices A and
    // B; the origin of the "G.2" label is unknown (not verified against any source), so it is
    // kept here only as a pointer to the old notes:
    //   tau(xi) ~ (L/Ledd) / Xi_tilde * (xi / h_d(xi))^2                         ("G.2")
    // Because n_e,b ~ xi'^-2 the integral is dominated by its lower limit, so tau decreases with
    // R where the disc flares (H/R grows -> the ray enters the corona farther out, where it is
    // less dense). z_disk = Height() the disc semi-height in the hot zone (up to Rhot) and
    // 0.01 R in the cold zone beyond Rhot (z_disk_to_R_cold_corona). Xi_tilde = 0.65
    // (T_ic/1e8 K)^-3/2 (Eq. 29b). Assumptions: isochoric corona (density n_e,b(xi) up to one
    // scale height, no vertical stratification of Eq. 34), corona hot (T_mb = T_ic), optically
    // thin, f = exp(-tau) + f_sc = 1 (attenuation and scattered radiation ignored); values
    // tau >~ 1 are outside the range of validity (the e^{-tau} factor in Eq. 38 would then
    // limit tau). R_IC is the genuine GM mu m_p/(k_B T_ic) (R_IC_corona(), no Dubus correction). Zero inside R_hc (OMK Eq. 20, see R_hc(): no hot corona there), at the
    // nodes shadowed by the disc itself (Shadow(), evaluated independently of Rhot: the cold disc
    // is typically in the shadow; the ray then does not reach the disc surface), and for wind
    // types without a T_ic.
    // NOTE: only the Compton-heated corona is included. Other outflows
    // (super-Eddington winds, magnetic winds) are ignored.
    if (!opt_str_.tau_r_corona) {
        vecd x(Nx(), 0.0);
        const double T_ic = wind_ ? wind_->T_ic() : 0.0;
        const double R_IC = R_IC_corona();
        if (T_ic > 0.0 && R_IC > 0.0) {
            const double R_hc_ = R_hc();
            const double l = Lbol_disk() / L_edd_disk();
            const double Xi_tilde = 0.65 * std::pow(T_ic / 1e8, -1.5);
            const vecd& shadow = Shadow();
            for (size_t i = first(); i < Nx(); i++) {
                const double z_disk = z_disk_corona(i);
                if (R()[i] < R_hc_ || shadow[i] > 0.0 || z_disk <= 0.0) {
                    continue;
                }
                // OMK91 Eq. (39) with the exact h_x = 0 solution of Eq. (40), xi_min = (z_disk/R)^2:
                x[i] = l / Xi_tilde * std::max(0.0, m::pow<2>(R()[i] / z_disk) - R_IC / R()[i]);
            }
        }
        opt_str_.tau_r_corona = std::move(x);
    }
    return *opt_str_.tau_r_corona;
}

double FreddiState::tau_r_corona_at_R_IC() const {
    // tau_r_corona (OMK91 Eqs. 39-40, exact for h_x = 0, see tau_r_corona()) evaluated directly
    // at the genuine Compton radius R_IC (R_IC_corona()), i.e. xi = 1, which is often outside
    // the radial grid: tau = (L/Ledd)/Xi_tilde * max(0, (R_IC / z_disk(R_IC))^2 - 1),
    // with z_disk = Height() (linearly interpolated) if R_IC <= Rhot
    // and z_disk_to_R_cold_corona * R_IC beyond Rhot. Zero if R_IC < R_hc or without a T_ic.
    // Shadow() is NOT taken into account here (it is a formula evaluation, not a grid value).
    const double T_ic = wind_ ? wind_->T_ic() : 0.0;
    const double R_IC = R_IC_corona();
    if (T_ic <= 0.0 || R_IC <= 0.0 || R_IC < R_hc()) {
        return 0.0;
    }
    double z_disk = z_disk_to_R_cold_corona * R_IC;
    if (R_IC <= R()[last()]) {
        const auto it = std::lower_bound(R().begin() + first(), R().begin() + last() + 1, R_IC);
        const size_t j = std::max<size_t>(it - R().begin(), first() + 1);
        const double w = (R_IC - R()[j-1]) / (R()[j] - R()[j-1]);
        z_disk = (1.0 - w) * Height()[j-1] + w * Height()[j];
    }
    if (z_disk <= 0.0) {
        return 0.0;
    }
    const double l = Lbol_disk() / L_edd_disk();
    const double Xi_tilde = 0.65 * std::pow(T_ic / 1e8, -1.5);
    return l / Xi_tilde * std::max(0.0, m::pow<2>(R_IC / z_disk) - 1.0);
}

size_t FreddiState::i_Fmax() const {
    // Grid index of the maximum of the viscous torque F (Fvis) over the hot zone, first()..last()
    // (up to Rhot), i.e. where Mdot = dF/dh = 0 (the boundary between the inner accreting and the
    // outer spreading disc). If F grows monotonically outward over the hot zone, this is Rhot.
    return std::max_element(F().begin() + first(), F().begin() + last() + 1) - F().begin();
}

// Optical depths from the central source to the disc surface at R()[i], used for the direct-flux
// attenuation f_dir = exp(-tau) in tau_perp_corona(). The name "tau_parallel" refers to the
// ray being nearly parallel to the disc plane; it is not OMK91's tau_parallel (Eqs. 48-49).
double FreddiState::tau_parallel_corona(size_t i) const {
    return tau_r_corona()[i];
}

double FreddiState::tau_parallel_wind(size_t i) const {
    return tau_r_wind()[i];
}

const vecd& FreddiState::tau_perp_corona() const {
    // Ostriker, McKee & Klein (1991, ApJ 377, 593), isothermal Compton-heated corona:
    //   tau_perp(xi) ~ f (L/Ledd) / (Xi_tilde xi^{1/2}) (T_mb/T_ic)^{1/2}      (O47, O55)
    //   Xi_tilde ~ 0.65 T_ic8^{-3/2}, T_ic8 = T_ic / 1e8 K                      (O29b)
    // xi = R / R_IC; T_mb/T_ic depends weakly on tau_parallel and is taken to be 1.
    // f is the local attenuation factor, Eqs. (17)-(18): f = f_dir + f_sc + f_d. As a first
    // approximation only f_dir = exp(-tau) is kept (OMK91, text after Eq. 45), with tau the
    // optical depth from the source to the disc surface at R: tau = tau_r_corona + tau_r_wind
    // (computed with f = 1 to avoid circularity; called tau_parallel below for the ray being
    // nearly parallel to the disc, not OMK91's tau_parallel of Eqs. 48-49).
    // R_IC is the genuine GM mu m_p/(k_B T_ic) (R_IC_corona()), not the Dubus-corrected
    // wind_->R_IC(). Zero inside R_hc (OMK Eq. 20, see R_hc()), at the nodes shadowed by the disc
    // itself (Shadow(), see tau_r_corona()) and for wind types without an inverse Compton
    // temperature. Other outflows (super-Eddington winds, magnetic winds) are not
    // included in tau_parallel.
    // WARNING / thin-corona cut: also zero where xi <= (z_disk/R)^2, i.e. where the corona scale
    // height h_c = xi^{3/2} R_IC (isothermal at T_ic, h_c/R = xi^{1/2}) is not above the disc
    // semi-height z_disk (z_disk_corona()). OMK's tau_perp (Eqs. 46-47) is derived for a corona
    // thicker than the disc (h_c > h_d); below that the isochoric model has no coronal layer and
    // Eq. (47) does not apply. This is the same condition, xi_min = (z_disk/R)^2 >= xi, that makes
    // tau_r_corona zero (exact h_x = 0 solution of Eq. 40, see tau_r_corona()). Without this cut
    // tau_perp_corona would be unattenuated (f = 1, since tau_r_corona = 0 there) and spuriously
    // large in the inner disc (e.g. ~0.4 at 3e7 cm in a test run where the cut is at ~3e9 cm).
    // NOTE: the corona there is not necessarily absent, only not described by this model (a
    // corona thinner than the disc, or the vertically stratified tail of Eq. 34).
    if (!opt_str_.tau_perp_corona) {
        vecd x(Nx(), 0.0);
        const double T_ic = wind_ ? wind_->T_ic() : 0.0;
        if (T_ic > 0.0 && R_IC_corona() > 0.0) {
            const double R_IC = R_IC_corona();
            const double R_hc_ = R_hc();
            const double l = Lbol_disk() / L_edd_disk();
            const double Xi_tilde = 0.65 * std::pow(T_ic / 1e8, -1.5);
            const vecd& shadow = Shadow();
            for (size_t i = first(); i < Nx(); i++) {
                const double z_disk = z_disk_corona(i);
                if (R()[i] < R_hc_ || shadow[i] > 0.0 || z_disk <= 0.0) {
                    continue;
                }
                const double xi = R()[i] / R_IC;
                // Thin-corona cut (see WARNING above): corona scale height <= disc semi-height
                if (xi <= m::pow<2>(z_disk / R()[i])) {
                    continue;
                }
                const double f = std::exp(-(tau_parallel_corona(i) + tau_parallel_wind(i)));
                x[i] = f * l / (Xi_tilde * std::sqrt(xi));
            }
        }
        opt_str_.tau_perp_corona = std::move(x);
    }
    return *opt_str_.tau_perp_corona;
}

const vecd& FreddiState::Mdot_wind_running() const {
    if (!opt_str_.Mdot_wind_running) {
        vecd mdot_wind(Nx(), 0.0);

        // Cumulative integration (trapezoidal rule) of dMdot_wind_dh() over the whole disc (not
        // just the hot region, see Mdot_wind()): Mdot_wind_running(h) = integral from h_in to h.
        mdot_wind[first()] = 0.0; // Boundary condition at Rin
        for (size_t i = first() + 1; i < Nx(); ++i) {
            const double dh = h()[i] - h()[i-1];
            const double avg_source = 0.5 * (dMdot_wind_dh(i) + dMdot_wind_dh(i-1));
            mdot_wind[i] = mdot_wind[i-1] + avg_source * dh;
        }

        opt_str_.Mdot_wind_running = std::move(mdot_wind);
    }
    return *opt_str_.Mdot_wind_running;
}


// double FreddiState::R_IC() const { 
//     return wind_ ? wind_->R_IC() : 0.0; 
// }

FreddiState::BasicWind::BasicWind(const FreddiState &state) :
    A_(state.Nx(), 0.),
    B_(state.Nx(), 0.),
    C_(state.Nx(), 0.),
    Vwind_(state.v_wind()), // state.v_wind() returns the esc velocity vector
    // R_IC_ (Compton radius) is only meaningful for wind models with an inverse Compton
    // temperature (Shields1986, Woods1996AGN, Woods1996ShieldsApproxWind); those recompute it
    // themselves in update() (see T_ic_current()), tracking windT_ic_approach. Other wind types
    // don't have "T_ic" in windparams at all, hence the guard here.
    R_IC_(state.args().disk->windparams.count("T_ic")
          ? (state.GM() * state.args().disk->mu * GSL_CONST_CGSM_MASS_PROTON) /
                (GSL_CONST_CGSM_BOLTZMANN * state.args().disk->windparams.at("T_ic"))
          : 0.0)
{}


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
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double R_g = 2*state.GM()/(m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT));
    const double L_edd = state.L_edd_disk();
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
        T_ic_const(state.args().disk->windparams.at("T_ic")),
        T_ic(T_ic_const),
        Pow(state.args().disk->windparams.at("Pow")) {
    update(state);
}

void FreddiState::Shields1986Wind::update(const FreddiState& state) {
    BasicWind::update(state);
    T_ic = state.T_ic_current(T_ic_const);
    T_ic_ = T_ic;
    const auto disk = state.args().disk;
    //  1983ApJ...271...70B page 4
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    // should the angular distribution be a factor here???
    //  1983ApJ...271...70B page 3
    const double R_iC = state.R_iC_from_T_ic(T_ic);
    R_IC_ = R_iC;
    const double L_edd = state.L_edd_disk();
    const double L_crit = (1.0 / 8.0) * std::sqrt(GSL_CONST_CGSM_MASS_ELECTRON / (disk->mu * GSL_CONST_CGSM_MASS_PROTON)) * std::sqrt((GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT ) / (GSL_CONST_CGSM_BOLTZMANN * T_ic)) * L_edd;
    const double P_iC = L / (4.0 * M_PI * m::pow<2>(R_iC) * Xi_max * GSL_CONST_CGSM_SPEED_OF_LIGHT);
    const double el = L/L_crit;
    

    // Integrated over the whole disc (not just the hot region): a cold, quiescent outer disc can
    // still carry a Compton wind if irradiated (see Mdot_wind()).
    for (size_t i = state.first(); i < state.Nx(); ++i) {
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
        T_ic_const(state.args().disk->windparams.at("T_ic")),
        T_ic(T_ic_const) {
    update(state);
}

void FreddiState::Woods1996AGNWind::update(const FreddiState& state) {
    BasicWind::update(state);
    T_ic = state.T_ic_current(T_ic_const);
    T_ic_ = T_ic;
    const auto disk = state.args().disk;
    const double L = state.Mdot_in() * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * state.eta();
    const double L_edd = state.L_edd_disk();
    const double R_iC = state.R_iC_from_T_ic(T_ic);
    R_IC_ = R_iC;
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


    // Integrated over the whole disc (not just the hot region): a cold, quiescent outer disc can
    // still carry a Compton wind if irradiated (see Mdot_wind()).
    for (size_t i = state.first(); i < state.Nx(); ++i) {
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
        T_ic_const(state.args().disk->windparams.at("T_ic")),
        T_ic(T_ic_const),
        Pow(state.args().disk->windparams.at("Pow")),
        IrAngDis(state.args().disk->windparams.at("IrAngDis")) {
    update(state);
}

void FreddiState::Woods1996ShieldsApproxWind::update(const FreddiState& state) {
    BasicWind::update(state);
    T_ic = state.T_ic_current(T_ic_const);
    T_ic_ = T_ic;
    const auto disk = state.args().disk;

    const double L = state.Lbol_disk();
    const double R_iC = state.R_iC_from_T_ic(T_ic);
    R_IC_ = R_iC;
    //const double VeL = std::sqrt(state.GM()/R_iC) ;
    //const double C_iC = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ic)/( GSL_CONST_CGSM_MASS_PROTON));
    //const double m_ch0 = disk->Mdot0 / (M_PI * state.R().back()*state.R().back());
    const double L_edd = state.L_edd_disk();
    const double L_crit = (1.0 / 8.0) * std::sqrt(GSL_CONST_CGSM_MASS_ELECTRON / (disk->mu * GSL_CONST_CGSM_MASS_PROTON)) * std::sqrt((GSL_CONST_CGSM_MASS_ELECTRON * GSL_CONST_CGSM_SPEED_OF_LIGHT* GSL_CONST_CGSM_SPEED_OF_LIGHT ) / (GSL_CONST_CGSM_BOLTZMANN * T_ic)) * L_edd;
    double el = L/L_crit;
    
    // Integrated over the whole disc (not just the hot region): a cold, quiescent outer disc can
    // still carry a Compton wind if irradiated (see Mdot_wind()).
    for (size_t i = state.first(); i < state.Nx(); ++i) {
        if (state.R()[i] > state.R_wind_inner_launch_radius(R_iC)) {
	    if (IrAngDis) {
			// Take account of the central flux angular distribution:
			el *= state.angular_dist_disk(state.Height()[i] / state.R()[i]) ; 
			//  angular_dist_disk(state.Height()[i] / state.R()[i]) ; 
			// disk_irr_source_->angular_dist(mu)
	    }
            //  1986ApJ...306...90S page 2
            const double xi = state.R()[i] / R_iC;
            const double xi1 = R_iC / state.R()[i];
            const double T_ch = T_ic * std::pow(el, 2.0 / 3.0) * std::pow(xi, -2.0 / 3.0);
            const double C_ch = std::sqrt((GSL_CONST_CGSM_BOLTZMANN * T_ch) / (disk->mu * GSL_CONST_CGSM_MASS_PROTON));
            const double C0 = (4.0 * M_PI * m::pow<3>(state.h()[i])) / (m::pow<2>(state.GM()));
            const double Fr =
                    L / (4.0 * M_PI * m::pow<2>(state.R()[i]) * Xi_max * C_ch * GSL_CONST_CGSM_SPEED_OF_LIGHT);

            const double Fc = ((std::pow((1 + m::pow<2>(((0.125 * el + 0.00382) * xi1))), (1.0 / 6.0))) /
                               std::pow((1 + m::pow<-2>((m::pow<4>(el)* (1.0 + 262.0 * m::pow<2>(xi))))),
                                        (1.0 / 6.0)));
	    
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

double FreddiState::v_visc(double r, double z2r_at_r) const {
	//oprel().Height(R()[i], F()[i])
	return  args().basic->alpha/oprel().Pi3 *  m::pow<2>(z2r_at_r)  * std::sqrt(GM()/r); 
	//Height2R()[]
}

double FreddiState::R_vis_struct (double r, double z2r_at_r) const {
        // estimate the radius, to which viscous process can reach in time tau
        return  r - v_visc(r, z2r_at_r) * args().calc->tau;       
}


double FreddiState::v_cooling_front(double r, double sigma_at_r) const {
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

double FreddiState::R_cooling_front (double r, double sigma_at_r) const {
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



void  FreddiState::find_R_max_where_Qirr_works () {
    // determine radius inside Rhot, outside of  which irradiation is  significant:
    // call  R_max_where_Qirr_works at the end of step in truncateOuterRadius
    set_maxR_Qirr_no_role (0.);
    for (size_t i = first(); i < last(); i++) {
	if ( ( Tirr().at(i) / Tph_vis().at(i) > args().disk->Tirr2Tvishot ) &&                                                    (maxR_Qirr_no_role() == 0.) ) { 
	    set_maxR_Qirr_no_role (R()[i]); 
	}  
    }
}
//current_.



double FreddiState::Tirr_critical (double r, int ii) const {
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
    double Qirr_Qvis = m::pow<4>(Tirr().at(ii) / Tph_vis().at(ii))*radius_popravka ;
    
    // determine radius, inside which irradiation is not significant:
    //find_R_max_where_Qirr_works();
    /*OLD   if ( (Qirr_Qvis < m::pow<4>(args().disk->Tirr2Tvishot) ) && (current_.maxR_Qirr_no_role == 0.) ) {
	set_maxR_Qirr_no_role (r); 
    }*/  
    
    double Tcrit; // value to return
    
    if (args().disk->Tirr_crit_approach == "const") { 
	 // since we calculate critical level for irradiation temperature, we use the fact that Tirr ~ R(-1/2) (This is abs true only if Cirr = const)
	Tcrit = args().disk->Thot * sqrt(radius_popravka);
	
	if (args().disk->scatter_by_corona == "no_") {    // @XRPCALCApril24
	    // do not take into account Tirr if Qirr/Qvis< critical_value
	    // and override previous assignment; instant return
	    if (Qirr_Qvis < m::pow<4>(args().disk->Tirr2Tvishot)) {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_Y Qirr/Qvis<crit \n" << std::endl;}
		// like the critical temperature is always too big 
		return 1e10;
	    } else
		// ( Tirr )
	    {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_W Qirr/Qvis=\n" <<Qirr_Qvis <<"\n"<<  std::endl;} 
	    }
	}
	
    } else if ((args().disk->Tirr_crit_approach == "Tavleev")  || (args().disk->Tirr_crit_approach == "Hameury") ) { 
        
	if (args().disk->scatter_by_corona == "no_") { //   @XRPCALCApril24
	    // if there is no scattering, the Rfront is in the shadow
	    // do not take into account Tirr if Qirr/Qvis < critical_value:
	    if (Qirr_Qvis < m::pow<4>(args().disk->Tirr2Tvishot)) {
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_V Qirr/Qvis=\n" <<Qirr_Qvis <<"\n"<<  std::endl;} 
		// cooling wave will be going with irradiation having no effect; 
		// return a large value of Tcrit
		// like the critical temperature is always too big 
		return 1e10;
	    }
	}

	if (Qirr_Qvis > 1.) { 
	        // since we calculate critical level for irradiation temperature, we use the fact that Tirr ~ R(-1/2) (This is completely true only if Cirr = const)
		Tcrit = (9040. - 2216.* 1./Qirr_Qvis ) * sqrt(radius_popravka) ;
		if (args().disk->Tirr_crit_approach == "Hameury") {
		    Tcrit = (8640. - 2216.* 1./Qirr_Qvis ) * sqrt(radius_popravka) ;
		}
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_F Qirr/Qvis>1 with Tcrit(Rhot)=" << Tcrit << " ii="<<ii<<"\n" << std::endl;}

	} else {
	        Tcrit = (9040. - 2216.) * pow(radius_popravka,0.5) ;
		if (args().calc->verb_level > VERB_LEVEL_MESSAGES) {std::cout << "c_G Qirr/Qvis<1 with Tcrit(Rhot)=" << Tcrit << " ii="<<ii<<"\n" << std::endl;}

	}
	if (args().disk->Tirr_crit_approach == "Hameury") {
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

#undef VERB_LEVEL_MESSAGES


