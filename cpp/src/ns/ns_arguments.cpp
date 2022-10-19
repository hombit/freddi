#include <algorithm> // max

#include <ns/ns_arguments.hpp>

#include <boost/math/tools/roots.hpp>

#include <fstream>

constexpr const char NeutronStarArguments::default_nsprop[];
constexpr const char NeutronStarArguments::default_Rm_type[];
constexpr const double NeutronStarArguments::default_h2r_bozzo;
constexpr const double NeutronStarArguments::default_chi_oblique;
constexpr const double NeutronStarArguments::default_hotspotarea;
constexpr const double NeutronStarArguments::default_epsilonAlfven;
//constexpr const double NeutronStarArguments::default_gamma_GL;
constexpr const double NeutronStarArguments::default_inversebeta;
constexpr const double NeutronStarArguments::default_Rdead;
constexpr const char NeutronStarArguments::default_fptype[];
constexpr const char NeutronStarArguments::default_kappat_type[];
constexpr const double NeutronStarArguments::default_kappat_value;
constexpr const char NeutronStarArguments::default_ns_grav_redshift[];
constexpr const double NeutronStarArguments::default_Rx_dummy;
constexpr const double NeutronStarArguments::default_freqx_dummy;

double NeutronStarArguments::initializeFreqx(const std::string& nsprop) {
	if (nsprop == "dummy") {
		return default_freqx_dummy;
	}
	if (nsprop == "sibgatullinsunyaev2000" || nsprop == "sibsun2000") {
		throw std::runtime_error("freqx must be specified for nsprop=sibgatullinsunyaev2000");
	}
	if (nsprop == "newt") {
		throw std::runtime_error("freqx must be specified for nsprop=newt");
	}
	throw std::invalid_argument("Wrong nsprop value");
}

double NeutronStarArguments::initializeRx(const std::string& nsprop, const std::optional<double> freqx) {
	if (nsprop == "dummy") {
		return default_Rx_dummy;
	}
	if (nsprop == "sibgatullinsunyaev2000" || nsprop == "sibsun2000") {
		if (!freqx) {
			throw std::runtime_error("freqx must be specified for nsprop=sibgatullinsunyaev2000");
		}
		return SibgatullinSunyaev2000Geometry::radiusNS(*freqx);
	}
	if (nsprop == "newt") {
		if (!freqx) {
			throw std::runtime_error("freqx must be specified for nsprop=newt");
		}
		return SibgatullinSunyaev2000Geometry::radiusNS(*freqx);
	}

	throw std::invalid_argument("Wrong nsprop value");
}

double NeutronStarArguments::R_Alfven(double GM, double Mdot) const {
	return epsilonAlfven * R_Alfven_basic(GM, Mdot);
}

double NeutronStarArguments::R_Alfven_basic(double GM, double Mdot) const {
	return std::pow(m::pow<4>(mu_magn) / (GM * m::pow<2>(Mdot)), 1./7.);
}

double NeutronStarArguments::R_Magn_KR07(double GM, double alpha, double Mdot) const {
	const double RA = R_Alfven_basic(GM, Mdot);
    const double RC = std::cbrt(GM / m::pow<2>(2*M_PI * freqx));
    //const double alpha = 0.5; //поправить потом  это дело
    const double chi_oblique_rad = chi_oblique * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2r_bozzo;
    
    double guess = 0.9;
    std::uintmax_t maxit = 300;
    double left = 0.;
    double right = 3.;
    boost::math::tools::eps_tolerance<double> tol(10);
    
    
    std::pair<double, double> r = boost::math::tools::toms748_solve(
            [RC, RA, chi_oblique_rad, parametr_bozzo, H_to_R_temp](double omega) 
			{ 
				return parametr_bozzo * (m::pow<2>(std::cos(chi_oblique_rad)) * (1. - omega) + 
				H_to_R_temp * (11. - 8. * omega) * m::pow<2>(std::sin(chi_oblique_rad)) ) * std::pow(RA/RC, 3.5) 
				- 0.5 * std::pow(omega, 10./3.); 
			},
                 left, right, tol, maxit);
				 std::cout << std::pow(r.first, 2./3.) * RC << std::endl;
    return std::pow(r.first, 2./3.) * RC;
    
}

double NeutronStarArguments::R_max_Fmagn_KR07(double GM, double alpha, double Mdot) const {
	const double RA = R_Alfven_basic(GM, Mdot);
    const double RC = std::cbrt(GM / m::pow<2>(2*M_PI * freqx));
    //const double alpha = 0.5; //поправить потом  это дело
    const double chi_oblique_rad = chi_oblique * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2r_bozzo;

	return std::pow((m::pow<2>(std::sin(chi_oblique_rad))*(11. * H_to_R_temp - 1. ) + 1.) 
	/ 
	(m::pow<2>(std::sin(chi_oblique_rad))*(8. * H_to_R_temp - 1. ) + 1.), 2./3.) * RC;
}

double NeutronStarArguments::R_Mdot_slope_KR07(double GM, double alpha, double Mdot) const {
	const double RA = R_Alfven_basic(GM, Mdot);
    const double RC = std::cbrt(GM / m::pow<2>(2*M_PI * freqx));
    //const double alpha = 0.5; //поправить потом  это дело
    const double chi_oblique_rad = chi_oblique * M_PI / 180.;
    const double parametr_bozzo = m::pow<2>(0.2) / alpha;
    const double H_to_R_temp = h2r_bozzo;
	const double Rmax = R_max_Fmagn_KR07(GM, alpha, Mdot);
	const double R0 = R_Magn_KR07(GM, alpha, Mdot); // We assume that Rmax > R0

	
    std::uintmax_t maxit = 300;
    double left = std::pow(RC/Rmax, 0.5);
    double right = std::pow(RC/R0, 0.5);
	double guess = left + 0.01;
    const int digits = std::numeric_limits<double>::digits/2;

	double result = boost::math::tools::newton_raphson_iterate(
            [RC, RA, chi_oblique_rad, parametr_bozzo, H_to_R_temp](double y) { 
				// y = 1/std::pow(R/RC, 0.5)
				double f = parametr_bozzo * (m::pow<2>(std::cos(chi_oblique_rad)) * (2*m::pow<10>(y) - 2*m::pow<7>(y))
				 + H_to_R_temp * (22*m::pow<10>(y) - 16.*m::pow<7>(y)) * m::pow<2>(std::sin(chi_oblique_rad)) ) * std::pow(RA/RC, 3.5) - 1;

				double f_deriv = parametr_bozzo * (m::pow<2>(std::cos(chi_oblique_rad)) * (20*m::pow<9>(y) - 14*m::pow<6>(y))
				 + H_to_R_temp * (220*m::pow<9>(y) - 16. * 7.*m::pow<6>(y)) * m::pow<2>(std::sin(chi_oblique_rad)) ) * std::pow(RA/RC, 3.5);

				return std::make_tuple(f, f_deriv); 
				},
                 guess, left, right, digits, maxit);
	return RC/m::pow<2>(result);
}
// написать все же через трансцендентное уравнение


double SibgatullinSunyaev2000Geometry::radiusNS(double freqx) {
	const double freqx_kHz = freqx / 1000.0;
	// Eq. 25
	const double r_km = 12.44 - 3.061 * freqx_kHz + 0.843 * m::pow<2>(freqx_kHz) + 0.6 * m::pow<3>(freqx_kHz)
	        + 1.56 * m::pow<4>(freqx_kHz);
	return 1e5 * r_km;
}

double SibgatullinSunyaev2000Geometry::radiusISCO(double freqx) {
	const double freqx_kHz = freqx / 1000.0;
	// Eq. 3, Eq. 26
	const double isco_minus_ns_km = 1.44 - 3.061 * freqx_kHz + 0.843 * m::pow<2>(freqx_kHz) + 0.6 * m::pow<3>(freqx_kHz) -
								 	0.22 * m::pow<4>(freqx_kHz);
	return 1e5 * isco_minus_ns_km + radiusNS(freqx);
}


NeutronStarBasicDiskBinaryArguments::NeutronStarBasicDiskBinaryArguments(
		const NeutronStarArguments& ns_args,
		double alpha_, std::optional<double> alphacold_,
		double Mx_, double kerr_,
		double period_,
		double Mopt_, double roche_lobe_fill_, double Topt_,
		std::optional<double> rin_, std::optional<double> rout_, std::optional<double> risco_
):
		BasicDiskBinaryArguments(
				alpha_, alphacold_,
				Mx_, kerr_,
				period_,
				Mopt_,
				roche_lobe_fill_,
				Topt_,
				rin_,
				rout_,
				initializeRisco(ns_args, risco_)) {
	if (!rin_) {
		rin = std::max(ns_args.Rx, risco);
	}
}

std::optional<double> NeutronStarBasicDiskBinaryArguments::initializeRisco(const NeutronStarArguments& ns_args, std::optional<double> risco){
	if (risco) {
		return risco;
	}
	if (ns_args.nsprop == "dummy") {
		return risco;
	}
	if (ns_args.nsprop == "newt") {
		return SibgatullinSunyaev2000Geometry::radiusISCO(ns_args.freqx);
        	}
	if (ns_args.nsprop == "sibgatullinsunyaev2000" || ns_args.nsprop == "sibsun2000") {
		return SibgatullinSunyaev2000Geometry::radiusISCO(ns_args.freqx);
	}
	throw std::invalid_argument("Wrong nsprop value");
}


NeutronStarDiskStructureArguments::NeutronStarDiskStructureArguments(
		const NeutronStarArguments& ns_args,
		const BasicDiskBinaryArguments &bdb_args,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot, double Tirr2Tvishot,
		const std::string& initialcond,
		std::optional<double> F0,
		std::optional<double> Mdisk0, std::optional<double> Mdot0,
		std::optional<double> powerorder,
		std::optional<double> gaussmu, std::optional<double> gausssigma,
		const std::string& wind, const pard& windparams
): DiskStructureArguments(opacity, OpacityRelated(opacity, bdb_args.Mx, bdb_args.alpha, mu),
						  Mdotout,
						  boundcond, Thot, Tirr2Tvishot,
						  initialcond,
						  initializeInitialFFunctionNS(OpacityRelated(opacity, bdb_args.Mx, bdb_args.alpha, mu),
													   ns_args, bdb_args,
													   initialcond, F0, Mdisk0, Mdot0,
													   powerorder,
													   gaussmu, gausssigma),
						  wind, windparams) {}


std::shared_ptr<DiskStructureArguments::InitialFFunction> NeutronStarDiskStructureArguments::initializeInitialFFunctionNS(
		const OpacityRelated& oprel,
		const NeutronStarArguments& ns_args,
		const BasicDiskBinaryArguments &bdb_args,
		const std::string& initialcond,
		std::optional<double> F0, std::optional<double> Mdisk0, std::optional<double> Mdot0,
		std::optional<double> powerorder,
		std::optional<double> gaussmu, std::optional<double> gausssigma) {
	if (!F0 && !Mdisk0 && !Mdot0) {
		throw std::runtime_error("One of F0, Mdisk0 or Mdot0 must be specified");
	}

	if (initialcond == "quasistat-ns" || initialcond == "quasistat_ns") {
		if (Mdot0) {
            double Ralfven;
            if (ns_args.Rm_type == "kluzniak" || ns_args.Rm_type == "Kluzniak"){
				double R0 = ns_args.R_Magn_KR07(bdb_args.Mx * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT, bdb_args.alpha, *Mdot0);
				/*double Rmax = ns_args.R_max_Fmagn_KR07(bdb_args.Mx * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT, *Mdot0);
				if (R0 >= Rmax) {
					Ralfven = R0;
				} else {
					Ralfven = std::max(ns_args.R_Mdot_slope_KR07(bdb_args.Mx * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT, *Mdot0), R0);
				}*/
				Ralfven = R0;
            } else {
				Ralfven = ns_args.R_Alfven( bdb_args.Mx * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT, *Mdot0);    
            }
			const double h_in = bdb_args.h(std::max(Ralfven, bdb_args.rin));
			const double h_out = bdb_args.h(bdb_args.rout);

			const double coeff = DiskStructureArguments::InitialFQuasistat::Coeff(h_in, h_out, oprel);
			F0 = *Mdot0 * (h_out - h_in) / h_out * h_in / oprel.f_F(h_in / h_out);
			Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;

			return std::make_shared<InitialFQuasistatNS>(*F0, *Mdisk0, *Mdot0, oprel, h_in);
		} else if (Mdisk0) {
			throw std::runtime_error("Mdisk0 is not supported for initialcond=quasistat-ns, please let us know if you need it");
		} else if (F0) {
			throw std::runtime_error("F0 is not supported for initialcond=quasistat-ns, please let us know if you need it");
		} else {
			throw std::logic_error("We couldn't be here");
		}
	}

	return DiskStructureArguments::initializeInitialFFunction(oprel, bdb_args, initialcond, F0, Mdisk0, Mdot0,
															  powerorder, gaussmu, gausssigma);
}


vecd NeutronStarDiskStructureArguments::InitialFQuasistatNS::operator()(const vecd& h) const {
	vecd F(h.size());
	for (size_t i = first(h); i < h.size(); ++i) {
		const double xi_LS2000 = h[i] / h.back();
		F[i] = F0 * oprel.f_F(xi_LS2000) * (1. - h_in / h[i]) / (1. - h_in / h.back());
		//F is Fsum, modification of F will be done in constructor of F. 
		//h[first(h)] does not equal h_in since h_in is calculated analitically, but h[first(h)] is the nearest mesh node.
		//It leads to the fact that F_vis(h[first(h)])!=0 in initial distribution.
	}
	return F;
}

size_t NeutronStarDiskStructureArguments::InitialFQuasistatNS::first(const vecd& h) const {
	size_t i;
	for (i = 0; h.at(i) < h_in; ++i) {}
	return i;
}

constexpr const char NeutronStarSelfIrradiationArguments::default_angular_dist_ns[];
