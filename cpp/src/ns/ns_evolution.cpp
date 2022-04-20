#include <exception>

#include <boost/math/tools/roots.hpp>

#include "exceptions.hpp"
#include "ns/ns_evolution.hpp"

#include <fstream>


FreddiNeutronStarEvolution::ConstKappaT::ConstKappaT(double value):
		value(value) {}

double FreddiNeutronStarEvolution::ConstKappaT::operator()(const FreddiNeutronStarEvolution& freddi, double R) const {
	return value;
}

FreddiNeutronStarEvolution::CorotationStepKappaT::CorotationStepKappaT(double in, double out):
		in(in), out(out) {}

double FreddiNeutronStarEvolution::CorotationStepKappaT::operator()(const FreddiNeutronStarEvolution& freddi, double R) const {
	if (R > freddi.R_cor()) {
		return out;
	}
	return in;
}

FreddiNeutronStarEvolution::Romanova2018KappaT::Romanova2018KappaT(double in, double out, double par1, double par2):
		in(in), out(out), par1(par1), par2(par2) {}

double FreddiNeutronStarEvolution::Romanova2018KappaT::operator()(const FreddiNeutronStarEvolution& freddi, double R) const {
	const double R_to_Rcor = R / freddi.R_cor();
	if (R_to_Rcor > 1.0) {
		// see Table 2 of Romanova+18 NewA 62, 94, the last line
		const double fastness = std::pow(R_to_Rcor, 1.5);
		const double feff = par1 * std::pow(fastness, par2);
		const double epsilonAlfven = freddi.epsilon_Alfven();
		const double out_minus_wind = out - feff * std::pow(epsilonAlfven, 3.5);
		if (out_minus_wind < 0){
			return 0.0;
		} else {
			return out_minus_wind;
		}
	}
	return in;
}


FreddiNeutronStarEvolution::NeutronStarStructure::NeutronStarStructure(
		const NeutronStarArguments &args_ns, FreddiEvolution* evolution):
		args_ns(args_ns),
		kappa_t(initialize_kappa_t(args_ns)),
		R_x(args_ns.Rx),
		redshift(initialize_redshift(evolution, args_ns)),
		R_m_min(std::max(R_x, evolution->R()[evolution->first()])),
		mu_magn(args_ns.mu_magn),
		R_cor(std::cbrt(evolution->GM() / m::pow<2>(2*M_PI * args_ns.freqx))),
		R_dead(args_ns.Rdead > 0. ? args_ns.Rdead : INFINITY),
//		F_dead((*kappa_t)(R_dead / R_cor) * m::pow<2>(mu_magn) / m::pow<3>(R_dead)),
		inverse_beta(args_ns.inversebeta),
		epsilon_Alfven(args_ns.epsilonAlfven),
		Rmtype(args_ns.Rm_type),
		h2rbozzo(args_ns.h2r_bozzo),
        chioblique(args_ns.chi_oblique),
		hot_spot_area(args_ns.hotspotarea),
		Fmagn(initialize_Fmagn(evolution)),
		dFmagn_dh(initialize_dFmagn_dh(evolution)),
		d2Fmagn_dh2(initialize_d2Fmagn_dh2(evolution)) {
	if (args_ns.Rdead > 0. && args_ns.Rdead < R_cor) {
		throw std::logic_error("R_dead is positive and less than R_cor, it is unacceptably");
	}
}


double FreddiNeutronStarEvolution::NeutronStarStructure::initialize_redshift(const FreddiEvolution* evolution, const NeutronStarArguments& args_ns) {
	if (args_ns.ns_grav_redshift == "off") {
		return 1.0;
	}
	if (args_ns.ns_grav_redshift == "on") {
		return 1.0 - 2.0 * evolution->R_g() / args_ns.Rx;
	}
	throw std::invalid_argument("Wrong nsgravredshift");
}


std::shared_ptr<FreddiNeutronStarEvolution::BasicKappaT> FreddiNeutronStarEvolution::NeutronStarStructure::initialize_kappa_t(const NeutronStarArguments& args_ns) {
	const auto& type = args_ns.kappat_type;
	const auto& params = args_ns.kappat_params;

	if (type == "const") {
		return std::make_shared<ConstKappaT>(params.at("value"));
	} else if (type == "corstep") {
		return std::make_shared<CorotationStepKappaT>(params.at("in"), params.at("out"));
	} else if (type == "romanova2018") {
		return std::make_shared<Romanova2018KappaT>(params.at("in"), params.at("out"), params.at("par1"), params.at("par2"));
	}
	throw std::invalid_argument("Wrong kappattype");
}


vecd FreddiNeutronStarEvolution::NeutronStarStructure::initialize_Fmagn(FreddiEvolution* evolution) const {
	vecd Fmagn_(evolution->Nx());
    
    if (Rmtype == "kluzniak" || Rmtype == "Kluzniak") {
        
		for (size_t i = 0; i < evolution->Nx(); i++) {
			Fmagn_[i] = F_Magn_KR07(evolution->R()[i]);
		}
	return Fmagn_;
    }
    
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
    
   if (Rmtype == "kluzniak" || Rmtype == "Kluzniak") {
        const double chi_oblique_rad = chioblique * M_PI / 180.;
        const double H_to_R_temp = h2rbozzo;
        //const double R_0 = std::max(R_Magn_KR07(evolution), R_x);
        const double k = inverse_beta * m::pow<2>(mu_magn) / std::pow(evolution->GM(), 0.5);
	double brackets1, brackets2;
	for (size_t i = 0; i < evolution->Nx(); i++) {
		brackets1 = - 2. / std::pow(evolution->R()[i], 3.5) + 2. * std::pow(R_cor, 1.5) / m::pow<5>(evolution->R()[i]);
        brackets2 = - 16./ (std::pow(evolution->R()[i], 3.5)) + 22. * std::pow(R_cor, 1.5) / m::pow<5>(evolution->R()[i]);
		dFmagn_dh_[i] = k * (m::pow<2>(std::cos(chi_oblique_rad)) * brackets1 + H_to_R_temp * m::pow<2>(std::sin(chi_oblique_rad)) * brackets2);
	}
	return dFmagn_dh_;
    }
    
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
    
    if (Rmtype == "kluzniak" || Rmtype == "Kluzniak") {
        const double chi_oblique_rad = chioblique * M_PI / 180.;
        const double H_to_R_temp = h2rbozzo;
        //const double R_0 = std::max(R_Magn_KR07(), R_x); здесь не нужно
        const double k = inverse_beta * m::pow<2>(mu_magn) / (evolution->GM());
	double brackets1, brackets2;
	for (size_t i = 0; i < evolution->Nx(); i++) {
		brackets1 = 14. /m::pow<4>(evolution->R()[i]) - 20. * std::pow(R_cor, 1.5) / std::pow(evolution->R()[i], 5.5);
        brackets2 = 112./m::pow<4>(evolution->R()[i]) - 220. * std::pow(R_cor, 1.5) / std::pow(evolution->R()[i], 5.5);
		d2Fmagn_dh2_[i] = k  * (m::pow<2>(std::cos(chi_oblique_rad)) * brackets1 + H_to_R_temp * m::pow<2>(std::sin(chi_oblique_rad)) * brackets2);
	}
	return d2Fmagn_dh2_;
    }
    
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

double FreddiNeutronStarEvolution::BasicNSMdotFraction::operator()(const FreddiNeutronStarEvolution& freddi, double R) const {
	const double Rx = freddi.R_x();
	const double Risco = freddi.args().basic->risco;

	if (R <= Rx || R <= Risco) {
		return 1.0;
	}

	const double Rcor = freddi.R_cor();
	return fp(R / Rcor);
}


double FreddiNeutronStarEvolution::NoOutflowNSMdotFraction::fp(double R_to_Rcor) const {
	return 1.;
}


double FreddiNeutronStarEvolution::PropellerNSMdotFraction::fp(double R_to_Rcor) const {
	return 0.;
}


double FreddiNeutronStarEvolution::CorotationBlockNSMdotFraction::fp(double R_to_Rcor) const {
	if (R_to_Rcor > 1.) {
		return 0.;
	}
	return 1.;
}


// https://arxiv.org/pdf/1010.1528.pdf Eksi-Kutlu (2010)
// fastness = (R_in/R_cor)^(3/2)  */
double FreddiNeutronStarEvolution::EksiKultu2010NSMdotFraction::fp(double R_to_Rcor) const {
	const double fastness2 = m::pow<3>(R_to_Rcor);
	double p = (1. - 1./fastness2);
	if (p < 0){
		p = 0;
	}
	return 1. - 1.5 * std::sqrt(p) + 0.5 * std::pow(p, 1.5);
}


FreddiNeutronStarEvolution::Romanova2018NSMdotFraction::Romanova2018NSMdotFraction(double par1, double par2):
		par1(par1), par2(par2) {}

// Return fraction of the accretion rate penetrating to the NS surface
// according to numerical results of MHD simulations (Romanova+2018, Table 2)
double FreddiNeutronStarEvolution::Romanova2018NSMdotFraction::fp(double R_to_Rcor) const {
	const double fastness = std::pow(R_to_Rcor, 1.5);
	double _fp = 0;
	if (fastness >= 1) {
		_fp = par1 * std::pow(fastness, par2);
	}
	if (_fp > 1) {
		_fp = 1.;
	}
	return 1. - _fp;
}


FreddiNeutronStarEvolution::GeometricalNSMdotFraction::GeometricalNSMdotFraction(double chi):
		chi(chi) {
	if (chi == 0) {
		throw std::logic_error("chi cannot be zero here, CorotationBlockNSMdotFraction should be used in this case");
	}
}

// Return fraction of the accretion rate penetrating to the NS surface
// for an inclined dipole with angle chi [grad] and R_to_R_cor = Rin/Rcor > 1
// see magnetospheric_form_1.mw
double FreddiNeutronStarEvolution::GeometricalNSMdotFraction::fp(double R_to_Rcor) const {
	const double mdot_factor = std::pow(R_to_Rcor, -3.5);
	if ( mdot_factor > 1. ) {
		return 1;
	}
	const double tmp = mdot_factor - m::pow<2>(std::cos(chi));
	if ( tmp < 0. ) {
		return 0;
	}
	return 1. - 2./M_PI * std::acos(std::sqrt(tmp) / std::sin(chi));
}


FreddiNeutronStarEvolution::BasicNSAccretionEfficiency::~BasicNSAccretionEfficiency() {}

double FreddiNeutronStarEvolution::BasicNSAccretionEfficiency::operator()(const FreddiNeutronStarEvolution& freddi, const double Rm) const {
	const double Rx = freddi.R_x();
	const double Risco = freddi.args().basic->risco;

	if ((Rx >= Risco) && (Rx >= Rm)) {
		return RxIsFurthest(freddi, Rm);
	}
	if ((Risco >= Rx) && (Risco >= Rm)) {
		return RiscoIsFurthest(freddi, Rm);
	}
	// ((Rm >= Rx) && (Rm >= Risco))
	// Rm could be NaN
	return RmIsFurthest(freddi, Rm);
}

double FreddiNeutronStarEvolution::DummyNSAccretionEfficiency::newtonian(const FreddiNeutronStarEvolution& freddi, const double Rm) const {
	const double R_in = freddi.R()[freddi.first()];
	const double Rg = freddi.R_g();
	const double Rx = freddi.R_x();
	return Rg * (1. / Rx - 0.5 / R_in);
}

double FreddiNeutronStarEvolution::RotatingNewtonianNSAccretionEfficiency::rotating_magnetosphere_newt(const FreddiNeutronStarEvolution& freddi, const double Rm)
const {
    const double Rsch = 2.0 * freddi.R_g();
    const double Rx = freddi.R_x();
    const double Rcor = freddi.R_cor();
    const double omega_ns = freddi.ns_str_->args_ns.freqx * 2 * M_PI;
    
    double Reff = Rm;
    // if the inner disc radius > Rcor, return efficiency calculated at R=Rcor:
    if (Rm > Rcor) {
    	Reff = Rcor;
    }
    const double omega_Kepl_Rin = std::sqrt(freddi.GM() / m::pow<3>(Reff));
    
    return Rsch / 2.0 / Rx * (1.0 - Rx / Reff)
    	+ m::pow<2>(omega_ns) / 2.0 / m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT) * (Rx-Reff) * (Rx+Reff)
    	+ Rsch / 4.0 / Reff * m::pow<2>(1.0 - omega_ns/omega_Kepl_Rin);
}

double FreddiNeutronStarEvolution::RotatingNewtonianNSAccretionEfficiency::small_magnetosphere_newt(const FreddiNeutronStarEvolution& freddi, const double Rm) const {
	const double Rsch = 2.0 * freddi.R_g();
	const double Rx = freddi.R_x();
	const double omega_ns = freddi.ns_str_->args_ns.freqx * 2 * M_PI;
	const double omega_Kepl_Rx = std::sqrt(freddi.GM() / m::pow<3>(Rx)); // Kepler angular velocity at NS surface

	return Rsch / 4.0 / Rx * m::pow<2>(1.0 - omega_ns/omega_Kepl_Rx);
}  


double FreddiNeutronStarEvolution::SibgatullinSunyaev2000NSAccretionEfficiency::schwarzschild(const FreddiNeutronStarEvolution& freddi, const double Rm) const {
	const double Rsch = 2.0 * freddi.R_g();
	const double Rx = freddi.R_x();
	return std::sqrt(1.0 - Rsch / Rm) - std::sqrt(1.0 - Rsch / Rx);
}

double FreddiNeutronStarEvolution::SibgatullinSunyaev2000NSAccretionEfficiency::rotating_magnetosphere_sibsun(const FreddiNeutronStarEvolution& freddi, const double Rm)
const {
    const double Rsch = 2.0 * freddi.R_g();
    const double Rx = freddi.R_x();
    const double Rcor = freddi.R_cor();
    const double omega_ns = freddi.ns_str_->args_ns.freqx * 2 * M_PI;
    const double omega_Kepl_Rx  = std::sqrt(freddi.GM() / m::pow<3>(Rx)); 
    const double etaSS2000 = small_magnetosphere(freddi, Rm);
    const double k = etaSS2000 / m::pow<2>(1.0 - omega_ns/omega_Kepl_Rx);
    // k =0.014 for Aql X-1
    
    double Reff = Rm;
    // if the inner disc radius > Rcor, return efficiency calculated at R=Rcor:
    if (Rm > Rcor) Reff = Rcor;
    const double omega_Kepl_Rin = std::sqrt(freddi.GM() / m::pow<3>(Reff)); 
    
    return std::sqrt(1.0 - Rsch / Reff) - std::sqrt(1.0 - Rsch / Rx)  
    + m::pow<2>(omega_ns) / 2.0 / m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT)  * (Rx-Reff) * (Rx+Reff)  
    + k * m::pow<2>(1.0 - omega_ns/omega_Kepl_Rin);
}

double FreddiNeutronStarEvolution::SibgatullinSunyaev2000NSAccretionEfficiency::small_magnetosphere(const FreddiNeutronStarEvolution& freddi, const double Rm) const {
	const double freqx_kHz = freddi.ns_str_->args_ns.freqx / 1000.0;
	// Eq. 1
	const double eta_ns_plus_disk = 0.213 - 0.153 * freqx_kHz + 0.02 * m::pow<2>(freqx_kHz);
	// Eq. 2
	const double ns_to_ns_plus_disk = 0.737 - 0.312 * freqx_kHz - 0.19 * m::pow<2>(freqx_kHz);
	const double eta_ns = eta_ns_plus_disk * ns_to_ns_plus_disk;
	if (eta_ns < 0.0 && eta_ns > 1.0) {
		throw std::logic_error("eta_ns must be in [0.0, 1.0]");
	}
	return eta_ns;
}


FreddiNeutronStarEvolution::FreddiNeutronStarEvolution(const FreddiNeutronStarArguments &args):
		FreddiEvolution(args),
		ns_str_(new NeutronStarStructure(*args.ns, this)),
		ns_irr_source_(initializeFreddiIrradiationSource(args.irr_ns->angular_dist_ns)),
		fp_(initializeNsMdotFraction(*args.ns)),
		eta_ns_(initializeNsAccretionEfficiency(*args.ns, this)) {
	// Change initial condition due presence of magnetic field torque. It can spoil user-defined initial disk
	// parameters, such as mass or Fout
       
	if (inverse_beta() <= 0.) {  // F_in is non-zero, Fmagn is zero everywhere
		current_.F_in = kappa_t(str_->R[0]) * m::pow<2>(mu_magn()) / m::pow<3>(R_cor());
		for (size_t i = 0; i < Nx(); i++) {
			current_.F[i] += current_.F_in;
		}
	} else {  // F_in is zero, and F + Fmagn = initial_cond + Fmagn_in
		//for (size_t i = 0; i < Nx(); i++) {             // before 13th December
            //current_.F[i] += -Fmagn()[i] + Fmagn()[0]; // before 13th December
          for (size_t i = first(); i < h().size(); ++i) {
            current_.F[i] += Fmagn()[first()] - Fmagn()[i];
//             std::cout << first() << std::endl;
//            std::cout << "Fmagn[" << i << "]=" << Fmagn()[i] << std::endl;
		}
	}
}

std::shared_ptr<FreddiNeutronStarEvolution::BasicNSMdotFraction> FreddiNeutronStarEvolution::initializeNsMdotFraction(const NeutronStarArguments& args_ns) {
	const auto& fptype = args_ns.fptype;
	const auto& fpparams = args_ns.fpparams;
	if (fptype == "no-outflow") {
		return std::make_shared<NoOutflowNSMdotFraction>();
	}
	if (fptype == "propeller") {
		return std::make_shared<PropellerNSMdotFraction>();
	}
	if (fptype == "corotation-block") {
		return std::make_shared<CorotationBlockNSMdotFraction>();
	}
	if (fptype == "eksi-kutlu2010") {
		return std::make_shared<EksiKultu2010NSMdotFraction>();
	}
	if (fptype == "romanova2018") {
		return std::make_shared<Romanova2018NSMdotFraction>(fpparams.at("par1"), fpparams.at("par2"));
	}
	if (fptype == "geometrical") {
		const double chi = fpparams.at("chi") * M_PI / 180.;
		if (chi == 0) {
			return std::make_shared<CorotationBlockNSMdotFraction>();
		} else {
			return std::make_shared<GeometricalNSMdotFraction>(chi);
		}
	}
	throw std::invalid_argument("Wrong fptype");
}


std::shared_ptr<FreddiNeutronStarEvolution::BasicNSAccretionEfficiency> FreddiNeutronStarEvolution::initializeNsAccretionEfficiency(const NeutronStarArguments& args_ns, const FreddiNeutronStarEvolution* freddi) {
	const auto& nsprop = args_ns.nsprop;
	if (nsprop == "dummy") {
		return std::make_shared<DummyNSAccretionEfficiency>();
	}
	if (nsprop == "sibgatullinsunyaev2000" || nsprop == "sibsun2000") {
		return std::make_shared<SibgatullinSunyaev2000NSAccretionEfficiency>();
	}
	if (nsprop == "newt") {
		return std::make_shared<RotatingNewtonianNSAccretionEfficiency>();
	}
	throw std::invalid_argument("Wrong nsprop");
}


double FreddiNeutronStarEvolution::Lbol_ns() const {
	return redshift() * Lbol_ns_rest_frame();
}


double FreddiNeutronStarEvolution::Lbol_ns_rest_frame() const {
        double Mdot = Mdot_in();
        if (Mdot < 0.0) Mdot = 0.0;
	return eta_ns() * fp() * Mdot * m::pow<2>(GSL_CONST_CGSM_SPEED_OF_LIGHT);
}


double FreddiNeutronStarEvolution::T_hot_spot() const {
	return std::pow(Lbol_ns_rest_frame() / (4*M_PI * hot_spot_area() * m::pow<2>(R_x()) * GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT), 0.25);
}


double FreddiNeutronStarEvolution::Lx_ns() {
	return redshift() * Lx_ns_rest_frame();
}


double FreddiNeutronStarEvolution::Lx_ns_rest_frame() {
	if (!ns_opt_str_.Lx_ns_rest_frame) {
		const double nu_min = args().flux->emin / redshift();
		const double nu_max = args().flux->emax / redshift();
		const double intensity = Spectrum::Planck_nu1_nu2(T_hot_spot(), nu_min, nu_max, 1e-4);
		ns_opt_str_.Lx_ns_rest_frame = 4*M_PI * hot_spot_area() * m::pow<2>(R_x()) * M_PI * intensity;
	}
	return *ns_opt_str_.Lx_ns_rest_frame;
}


void FreddiNeutronStarEvolution::invalidate_optional_structure() {
	FreddiEvolution::invalidate_optional_structure();
	ns_opt_str_ = NeutronStarOptionalStructure();
}


void FreddiNeutronStarEvolution::truncateInnerRadius() {
 
    //const auto& Rm_definition = args_ns.Rm_definition;
	if (R_dead() <= 0.) {
		return;
	}
	
	//std::cout << Fmagn()[current_.first] << " ***************** " << Mdot_in() << " " << Mdot_in_prev() << std::endl;
	// check proper value of accretion rate:)
	if (((!std::isfinite(Mdot_in_prev())) || ( Mdot_in() < 0.0 ) || ( Mdot_in_prev() < 0.0 ))) { 
		return;
	}
	// check that Mdot decaying:
	if (( Mdot_in() > Mdot_in_prev())  )  { //&& inverse_beta() <= 0.
		return;
	}
    double R_magnetic;
    if (Rmtype() == "bozzo" || Rmtype() == "Bozzo"){
        R_magnetic = R_Magn_bozzo18(); 
    }
    
    if (Rmtype() == "alfven" || Rmtype() == "Alfven"){
        R_magnetic = R_Alfven();
    }
    
    if (Rmtype() == "Kluzniak" || Rmtype() == "kluzniak"){
		//std::cout << R_Magn_KR07() << "*" << R_max_Fmagn_KR07() << std::endl;
        //R_magnetic = std::max(R_Magn_KR07(), R_Mdot_slope_KR07());
		R_magnetic = R_Magn_KR07();
    }

    double R_m = std::max(R_m_min(), R_magnetic);
    R_m = std::min(R_m, R_dead());
    /*
    if (R_m == R_magnetic){
        std::cout << "defined my magnetosphere\n";
    }
    
    if (R_m == R_m_min()){
        std::cout << "defined my NS\n";
    }
    
    if (R_m == R_dead()){
        std::cout << "defined my dead disk\n";
    }
    //std::cout << Mdot_in() << std::endl;
    */
    size_t ii;
    for (ii = first(); ii <= last()-2; ii++) {
        if (R()[ii + 1] > R_m){
            break;
        }
    }
    if (ii >= last() - 2) {
        throw RadiusCollapseException();
    }

	double R0_KR07 = R_magnetic;
	double Mdot_KR07 = Mdot_in();

	if (Rmtype() == "Kluzniak" || Rmtype() == "kluzniak"){
		if ((R()[ii+1] >= R_max_Fmagn_KR07()) || (ii == first())) {
    		R_m = R()[ii];
    		current_.first = ii;
		} else {
			R_m = R()[ii+1];
    		current_.first = ii+1;
		}
	} else {
		R_m = R()[ii];
    	current_.first = ii;
	}
    double new_F_in = 0;
        if (inverse_beta() <= 0.) {
            if (R_m <= R_cor()) {
                new_F_in = kappa_t(R_m) * m::pow<2>(mu_magn()) / m::pow<3>(R_cor());
            } else {
    //			new_F_in = F_dead() * m::pow<3>(R_dead() / R_m);
                new_F_in = kappa_t(R_m) * m::pow<2>(mu_magn()) / m::pow<3>(R_m);
            }
        } else {
			if (Rmtype() == "Kluzniak" || Rmtype() == "kluzniak"){
            	//new_F_in = Mdot_KR07*(std::pow(GM()*R_m, 0.5)-std::pow(GM()*R0_KR07, 0.5)) - F_Magn_KR07(R_m) + F_Magn_KR07(R0_KR07);
				/*std::cout << R0_KR07 << " " << R_m << " " << new_F_in << " " << Mdot_KR07 << " " 
				<< Mdot_KR07*(std::pow(GM()*R_m, 0.5)-std::pow(GM()*R0_KR07, 0.5)) <<  " " <<  F_Magn_KR07(R_m) - F_Magn_KR07(R0_KR07) << " " <<
				(Mdot_KR07*(std::pow(GM()*R_m, 0.5)-std::pow(GM()*R0_KR07, 0.5)))/(F_Magn_KR07(R_m) - F_Magn_KR07(R0_KR07)) << std::endl;
				*/
			new_F_in = 0;
			} else {
				new_F_in = 0;
			}
        }
        current_.F_in = new_F_in;
}
double FreddiNeutronStarEvolution::Mdot_in() const {
	const double dF_dh = (F()[first() + 1] - F()[first()]) / (h()[first() + 1] - h()[first()]);
	return dF_dh + dFmagn_dh()[first()]; //F is F_vis
    //return dF_dh;
}

double FreddiNeutronStarEvolution::R_Alfven() const {
	return ns_str_->args_ns.R_Alfven(GM(), Mdot_in());
}

double FreddiNeutronStarEvolution::R_Alfven_basic() const {
	return ns_str_->args_ns.R_Alfven_basic(GM(), Mdot_in());
}


double FreddiNeutronStarEvolution::R_Magn_bozzo18() const {
	const double chi_oblique_rad = chioblique() * M_PI / 180.;
    const double alpha = args().basic->alpha; 
	//eta (Bozzo2018) = 0.2, screening factor
	const double parametr_bozzo = 2. * m::pow<2>(0.2) / alpha; 
	const double RA = R_Alfven_basic();
    const double RCorrot = R_cor();
	const double parametr_thetta = RCorrot  / RA;


    const double H_to_R_temp = h2rbozzo(); //временно
	const bool sign_at_inf = 1; 
    
    //sign of Bozzo expression should be minus in infunuty!!!

    
    const double xbozzo_0 = R()[first()] / RCorrot;
    const double fxbozzo_0 = parametr_bozzo * ( (1. - std::pow(xbozzo_0, 1.5) ) * m::pow<2>(std::cos(chi_oblique_rad)) 
	+ H_to_R_temp * (8. - 5. * std::pow(xbozzo_0, 1.5)) * m::pow<2>(std::sin(chi_oblique_rad))) - std::pow(xbozzo_0 * parametr_thetta, 3.5);
    
    if (std::signbit(fxbozzo_0) == sign_at_inf){//that would mean the root is smaller than R_m_min and so we return any number less than NS radius 
        return 0. ; 
    }
    
	size_t jj;
    double xbozzo;
    double fxbozzo;
    //double fxbozzo_prev;
    //now compare the sign of expression with its sign at the infinity, break when sign changes
    
	for(jj = first(); jj <= last() - 2 ; jj++){
        //double fxbozzo_prev = fxbozzo;
		xbozzo = R()[jj] / RA;
		fxbozzo = parametr_bozzo *( (1. - std::pow(xbozzo / parametr_thetta, 1.5) ) * m::pow<2>(std::cos(chi_oblique_rad)) 
		+ H_to_R_temp * (8. - 5. * std::pow(xbozzo / parametr_thetta, 1.5)) * m::pow<2>(std::sin(chi_oblique_rad))) - std::pow(xbozzo, 3.5);
        if (fxbozzo < 0.){ //check if sigh of expression changes
            break;
        }
	}
	if (jj == last() - 2){
        throw RadiusCollapseException();;
    }

	return R()[jj];
	//мы могли бы провести здесь интерполяцию между R[jj_найденное] и R[jj_прошлое], тогда, наверное, даже не будет ёлочки. 
    //return (R[jj] * fxbozzo_prev - R[jj + 1] * fxbozzo) / (fxbozzo_prev - fxbozzo)
}

/*double FreddiNeutronStarEvolution::R_Magn_KR07() const {
	const double RA = R_Alfven_basic();
    const double RC = R_cor();
    const double alpha = 0.5; //поправить потом  это дело
    const double chi_oblique_rad = chioblique() * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2rbozzo();

	double guess = 0.9;
    std::uintmax_t maxit = 1000;
    double left = 0.;
    double right = 3.;
    boost::math::tools::eps_tolerance<double> tol(10);
    
    
    std::pair<double, double> r = boost::math::tools::toms748_solve(
            [RC, RA, chi_oblique_rad, parametr_bozzo, H_to_R_temp](double omega) 
			{ 
				return parametr_bozzo * (m::pow<2>(std::cos(chi_oblique_rad)) * (1. - omega) + 
			H_to_R_temp * (11. - 8. * omega) * m::pow<2>(std::sin(chi_oblique_rad)) ) * 
			std::pow(RA/RC, 3.5) - 0.5 * std::pow(omega, 10./3.);
			},
                 left, right, tol, maxit);
	double omega = r.first;
    return std::pow(r.first, 2./3.) * RC;
	//return ns_str_->args_ns.R_Magn_KR07(GM(), Mdot_in()); 
}*/

double FreddiNeutronStarEvolution::F_Magn_KR07(const double R) const {
	const double chi_oblique_rad = chioblique() * M_PI / 180.;
	const double H_to_R_temp = h2rbozzo();
//  const double R_0 = std::max(R_Magn_KR07(), R_x); не нужно сейчас
	const double k = inverse_beta() * m::pow<2>(mu_magn()) / 9.;
	double brackets1, brackets2;
	double Fmagn;
	brackets1 = 3. / m::pow<3>(R) - 2. * std::pow(R_cor(), 1.5) / std::pow(R, 4.5);
	brackets2 = 24. / m::pow<3>(R) - 22. * std::pow(R_cor(), 1.5) / std::pow(R, 4.5);
	Fmagn = k  * (m::pow<2>(std::cos(chi_oblique_rad)) * brackets1 + H_to_R_temp * m::pow<2>(std::sin(chi_oblique_rad)) * brackets2);
	return Fmagn;
    
}

/*double FreddiNeutronStarEvolution::NeutronStarStructure::R_Magn_KR07(FreddiEvolution* evolution) const {
	const double RA = R_Alfven_basic;
    const double RC = R_cor;
    const double alpha = 0.5; //поправить потом  это дело
    const double chi_oblique_rad = chioblique * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2rbozzo;

	double guess = 0.9;
    std::uintmax_t maxit = 300;
    double left = 0.;
    double right = 3.;
    boost::math::tools::eps_tolerance<double> tol(10);
    
    
    std::pair<double, double> r = boost::math::tools::toms748_solve(
            [RC, RA, chi_oblique_rad, parametr_bozzo, H_to_R_temp](double omega) 
			{ 
				return parametr_bozzo * (m::pow<2>(std::cos(chi_oblique_rad)) * (1 - omega) + 
			H_to_R_temp * (11 - 8 * omega) * m::pow<2>(std::sin(chi_oblique_rad)) ) * std::pow(RA/RC, 3.5) 
			- 0.5 * std::pow(omega, 10./3.); 
			},
                 left, right, tol, maxit);
    return std::pow(r.first, 2./3.) * RC;
	//return args_ns.R_Magn_KR07(evolution->GM(), evolution->Mdot_in());
    
}*/

double FreddiNeutronStarEvolution::R_max_Fmagn_KR07() const { // R : Fmagn -> max
    const double RC = R_cor(); 
    const double chi_oblique_rad = chioblique() * M_PI / 180.;
    const double H_to_R_temp = h2rbozzo();
    
    return std::pow((m::pow<2>(std::sin(chi_oblique_rad))*(11 * H_to_R_temp - 1 ) + 1) 
	/ 
	(m::pow<2>(std::sin(chi_oblique_rad))*(8.* H_to_R_temp - 1 ) + 1), 2./3.) * RC;
    
}

double FreddiNeutronStarEvolution::R_Mdot_slope_KR07() const { // R : dFmag_dh = Mdot 
    const double RA = R_Alfven_basic();                        // Not used now
    const double RC = R_cor();
	const double R_max = R_max_Fmagn_KR07();
    const double alpha = args().basic->alpha; 
    const double chi_oblique_rad = chioblique() * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2rbozzo();

	size_t ii;
    for (ii = first(); ii <= last() - 2; ii++) {
        if (R()[ii + 1] > R_max){
            break;
        }
    }
    if (ii >= last() - 2) {
        throw RadiusCollapseException();
    }

	size_t jj;
	//We assume that dFmagn_dh < Mdot

	for (jj = ii; jj >= first(); jj--) {
        if (dFmagn_dh()[jj] >= Mdot_in()) {
			break;
		}
    }

    return R()[jj];
    
}

double FreddiNeutronStarEvolution::NeutronStarStructure::F_Magn_KR07(const double R) const {
        const double chi_oblique_rad = chioblique * M_PI / 180.;
        const double H_to_R_temp = h2rbozzo;
//      const double R_0 = std::max(R_Magn_KR07(), R_x); не нужно сейчас
        const double k = inverse_beta * m::pow<2>(mu_magn) / 9.;
	    double brackets1, brackets2;
        double Fmagn;
        brackets1 = 3. / m::pow<3>(R) - 2. * std::pow(R_cor, 1.5) / std::pow(R, 4.5);
        brackets2 = 24. / m::pow<3>(R) - 22. * std::pow(R_cor, 1.5) / std::pow(R, 4.5);
		Fmagn = k  * (m::pow<2>(std::cos(chi_oblique_rad)) * brackets1 + H_to_R_temp * m::pow<2>(std::sin(chi_oblique_rad)) * brackets2);
        return Fmagn;
    
}


IrradiatedStar::sources_t FreddiNeutronStarEvolution::star_irr_sources() {
	auto sources = FreddiEvolution::star_irr_sources();
	sources.push_back(ns_irr_source_->irr_source(*this, Lbol_ns()));
	return sources;
}


vecd FreddiNeutronStarEvolution::windC() const {
	auto C = FreddiEvolution::windC();
	for (size_t i = 0; i < C.size(); i++) {
		C[i] += d2Fmagn_dh2()[i];
	}
	return C;
}


const vecd& FreddiNeutronStarEvolution::Qx() {
	if (!opt_str_.Qx) {
		vecd x(Nx());
		const vecd& K = Kirr();
		const vecd& H = Height();
		const double L_disk = Lbol_disk();
		const double L_ns = Lbol_ns();
		for (size_t i = first(); i < Nx(); i++) {
			const double mu = H[i] / R()[i];
			x[i] = K[i] * (L_disk * angular_dist_disk(mu) + L_ns * angular_dist_ns(mu)) / (4. * M_PI * m::pow<2>(R()[i]));
		}
		opt_str_.Qx = std::move(x);
	}
	return *opt_str_.Qx;
}

/*
const vecd& FreddiNeutronStarEvolution::Tph_X() {
	if (!opt_str_.Tph_X) {
		vecd x(Nx());
		const vecd& Tvis = Tph_vis();
		for (size_t i = first(); i < Nx(); i++) {
			x[i] = args().flux->colourfactor * Tvis[i];
		}
		opt_str_.Tph_X = std::move(x);
	}
	return *opt_str_.Tph_X;
}
*/

double FreddiNeutronStarEvolution::Lbol_disk() const {
	return (F()[first()] + 0.5 * Mdot_in() * h()[first()]) * omega_i(first());
}
