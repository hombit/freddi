#include <boost/numeric/odeint.hpp>

#include "arguments.hpp"
#include "util.hpp"


namespace odeint = boost::numeric::odeint;


constexpr const char GeneralArguments::default_prefix[];
constexpr const char GeneralArguments::default_dir[];
constexpr const unsigned short GeneralArguments::default_output_precision;
constexpr const unsigned int GeneralArguments::default_temp_sparsity_output;


constexpr const double BasicDiskBinaryArguments::default_alpha_to_alphacold;
constexpr const double BasicDiskBinaryArguments::default_kerr;
constexpr const double BasicDiskBinaryArguments::default_roche_lobe_fill;
constexpr const double BasicDiskBinaryArguments::default_Topt;


constexpr const char DiskStructureArguments::default_opacity[];
constexpr const double DiskStructureArguments::default_Mdotout;
constexpr const char DiskStructureArguments::default_initialcond[];
constexpr const double DiskStructureArguments::default_Thot;
constexpr const char DiskStructureArguments::default_boundcond[];
constexpr const char DiskStructureArguments::default_scatter_by_corona[];
constexpr const char DiskStructureArguments::default_DIM_front_approach[];
constexpr const double DiskStructureArguments::mu;
constexpr const char DiskStructureArguments::default_wind[];

DiskStructureArguments::DiskStructureArguments(
		const BasicDiskBinaryArguments &bdb_args,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot, double Tirr2Tvishot,
		double Rfront_Mdotzero_factor, double DIM_front_Mdot_factor, const std::string& check_state_approach, 
		const std::string& check_Sigma_approach,
		const std::string& check_Temp_approach,
		const std::string& DIM_front_approach,
		const std::string& scatter_by_corona,
		const std::string& initialcond,
		std::optional<double> F0,
		std::optional<double> Mdisk0, std::optional<double> Mdot0,
		std::optional<double> powerorder,
		std::optional<double> gaussmu, std::optional<double> gausssigma,
		const std::string& wind, const pard& windparams
	):
		opacity(opacity),
		oprel(OpacityRelated(opacity, bdb_args.Mx, bdb_args.alpha, mu)),
		Mdotout(Mdotout),
		boundcond(boundcond),
		Thot(Thot),
		Tirr2Tvishot(Tirr2Tvishot),
		Rfront_Mdotzero_factor(Rfront_Mdotzero_factor), DIM_front_Mdot_factor(DIM_front_Mdot_factor), check_state_approach(check_state_approach),check_Sigma_approach(check_Sigma_approach), check_Temp_approach(check_Temp_approach),
		DIM_front_approach(DIM_front_approach),scatter_by_corona(scatter_by_corona),
		initialcond(initialcond),
		initial_F_function(initializeInitialFFunction(oprel,
													  bdb_args, initialcond,
													  F0, Mdisk0, Mdot0,
													  powerorder,
													  gaussmu, gausssigma)),
		wind(wind),
		windparams(windparams) {}

std::shared_ptr<DiskStructureArguments::InitialFFunction> DiskStructureArguments::initializeInitialFFunction(
				const OpacityRelated& oprel,
				const BasicDiskBinaryArguments &bdb_args, const std::string& initialcond,
				std::optional<double> F0, std::optional<double> Mdisk0, std::optional<double> Mdot0,
				std::optional<double> powerorder,
				std::optional<double> gaussmu, std::optional<double> gausssigma) {
		if (!F0 && !Mdisk0 && !Mdot0) {
			throw std::runtime_error("One of F0, Mdisk0 or Mdot0 must be specified");
		}

		const double h_in = bdb_args.h(bdb_args.rin);
		const double h_out = bdb_args.h(bdb_args.rout);

		if (initialcond == "linearF" || initialcond == "linear") {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel.m;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			const double coeff = std::pow(h_out - h_in, oprel.n + 1.) * integral / ((1. - oprel.m) * oprel.D);

			if (Mdot0) {
				F0 = *Mdot0 * (h_out - h_in);
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			} else if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
				Mdot0 = *F0 / (h_out - h_in);
			} else if (F0) {
				Mdot0 = *F0 / (h_out - h_in);
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			} else {
				throw std::logic_error("We couldn't be here");
			}

			return std::make_shared<InitialFLinearF>(*F0, *Mdisk0, *Mdot0);
		}
		if (initialcond == "powerF" || initialcond == "power") {
			if (!powerorder) {
				throw std::runtime_error("powerorder must be specified for the case of initialcond=powerF");
			}

			if (Mdot0) {
				throw std::runtime_error("Set Mdisk0 or F0 instead of Mdot0 if initialcond=powerF");
			}

			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = (1. - oprel.m) * *powerorder;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			const double coeff = std::pow(h_out - h_in, oprel.n + 1.) * integral / ((1. - oprel.m) * oprel.D);

			if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
			} else {
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			}

			return std::make_shared<InitialFPowerF>(*F0, 0.0, *Mdisk0, *powerorder);
		}

		if (initialcond == "powerSigma") {
			if (!powerorder) {
				throw std::runtime_error("powerorder must be specified for the case of initialcond=powerSigma");
			}

			if (Mdot0) {
				throw std::runtime_error("Set Mdisk0 or F0 instead of Mdot0 if initialcond=powerSigma");
			}

			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = *powerorder;
			const double b = 3.;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			const double coeff = m::pow<4>(h_out - h_in) * integral /
								 ((1. - oprel.m) * oprel.D * std::pow(h_out, 3. - oprel.n));

			if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
			} else {
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			}

			return std::make_shared<InitialFPowerSigma>(*F0, 0.0, *Mdisk0, *powerorder, oprel);
		}

		if (initialcond == "sineF" || initialcond == "sine" ||
			initialcond == "sinusF" || initialcond == "sinus") {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel.m;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(sin(x * M_PI_2), a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			const double coeff = pow(h_out - h_in, oprel.n + 1.) * integral / ((1. - oprel.m) * oprel.D);

			if (Mdot0) {
				F0 = *Mdot0 * (h_out - h_in) * 2. / M_PI;
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			} else if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
				Mdot0 = *F0 * M_PI / (2. * (h_out - h_in));
			} else if (F0) {
				Mdot0 = *F0 * M_PI / (2. * (h_out - h_in));
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			} else {
				throw std::logic_error("We couldn't be here");
			}

			return std::make_shared<InitialFSineF>(*F0, *Mdisk0, *Mdot0);
		}

		if (initialcond == "gaussF") {
			if (!gaussmu) {
				throw std::runtime_error("gaussmu must be specified for the case of initialcond=gaussF");
			}
			if (!gausssigma) {
				throw std::runtime_error("gausssigma must be specified for the case of initialcond=gaussF");
			}

			if (*gaussmu <= 0. or *gaussmu > 1.) {
				throw std::runtime_error("gaussmu value should be in [0..1]");
			}

			if (Mdot0) {
				throw std::runtime_error(
						"Setting Mdot if initialcond=gaussF produces unstable results and it isn't motivated physically. Set F0 or Mdisk0 instead");
				//F0_ = Mdot_in * (h_out - h_in) * m::pow<2>(gauss_sigma) / gauss_mu * exp( m::pow<2>(gauss_mu) / (2. * m::pow<2>(gauss_sigma)) );
			}

			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel.m;
			const double b = oprel.n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0, gaussmu, gausssigma](const double &y, double &dydx, double x) {
						dydx = exp(-m::pow<2>(x - *gaussmu) * a / (2. * m::pow<2>(*gausssigma))) *
							   pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			const double coeff = std::pow(h_out - h_in, oprel.n + 1.) * integral / ((1. - oprel.m) * oprel.D);

			if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
			} else {
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			}

			return std::make_shared<InitialFGaussF>(*F0, 0.0, *Mdisk0, *gaussmu, *gausssigma);
		}

		if (initialcond == "quasistat") {
			const double coeff = InitialFQuasistat::Coeff(h_in, h_out, oprel);
			// How to find F0 from Mdot0?
			//
			// Let f1 be the actual F(h) distribution in the disc.
			// Then Mdot = dF/dh = F0 / hout * f1'(xi), xi = h/h_out,
			// and Mdot0 = Mdot (h=h_in) = F0 / hout * f1'(xi_in), where xi_in = h_in/h_out
			//
			// For f1 (xi) = f_F(xi) * (1-h_in/h_out/xi) / (1-h_int/h_out) 
			// (see Lipunova & Malanchev (2012), Appendix B, case (v))
			// it can be shown that f1'(xi_in) = f_F(xi_in) * h_out / h_in / (1-h_in/h_out)
			// Thus F0 = Mdot0 * h_in/h_out * (h_out-h_in) / f_F(xi_in)
			//
			if (Mdot0) {
				F0 = *Mdot0 * (h_out - h_in) / h_out * h_in / oprel.f_F(h_in / h_out);
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
			} else if (Mdisk0) {
				F0 = std::pow(*Mdisk0 / coeff, 1. / (1. - oprel.m));
				Mdot0 = *F0 * h_out / (h_out - h_in) * oprel.f_F(h_in / h_out) / h_in;
			} else if (F0) {
				Mdisk0 = std::pow(*F0, 1. - oprel.m) * coeff;
				Mdot0 = *F0 * h_out / (h_out - h_in) * oprel.f_F(h_in / h_out) / h_in;
			} else {
				throw std::logic_error("We couldn't be here");
			}

			return std::make_shared<InitialFQuasistat>(*F0, *Mdisk0, *Mdot0, oprel);
		}

		throw std::runtime_error("Invalid value of initialcond");
}

DiskStructureArguments::InitialFFunction::~InitialFFunction() {}

vecd DiskStructureArguments::InitialFPowerF::operator()(const vecd& h) const {
	vecd F(h.size());
	for (size_t i = 0; i < h.size(); ++i) {
		F[i] = F0 * std::pow((h[i] - h.front()) / (h.back() - h.front()), powerorder);
	}
	return F;
}

vecd DiskStructureArguments::InitialFPowerSigma::operator()(const vecd& h) const {
	vecd F(h.size());
	for (size_t i = 0; i < h.size(); ++i) {
		const double Sigma_to_Sigmaout = std::pow((h[i] - h.front()) / (h.back() - h.front()), powerorder);
		F[i] = F0 * std::pow(h[i] / h.back(), (3. - oprel.n) / (1. - oprel.m)) *
				std::pow(Sigma_to_Sigmaout, 1. / (1. - oprel.m));
	}
	return F;
}

vecd DiskStructureArguments::InitialFSineF::operator()(const vecd& h) const {
	vecd F(h.size());
	for (size_t i = 0; i < h.size(); ++i) {
		F[i] = F0 * std::sin((h[i] - h.front()) / (h.back() - h.front()) * M_PI_2);
	}
	return F;
}

vecd DiskStructureArguments::InitialFQuasistat::operator()(const vecd& h) const {
	// As an approximation for the initial distribution of F(h), the following function
	// is adopted : f1 (xi) = f_F(xi) * (1-h_in/h_out/xi) / (1-h_int/h_out), where 
	// xi = h/h_out and
	// f_F(xi) is the dimensionless function found by Lipunova & Shakura (2000)
	// as a part of the self-similar solution of the viscous diffusion equation 
	// in the interval [0..h_out]
	// see Appendix B of Lipunova & Malanchev (2017), case (v).
	// f1 is not a solution of the diffusion equation since there is no 
	// (or we did not find) self-similar solution in the case h_in > 0.
	// However, f1 serves as an initial distribution quite well; 
	// it satisfies f1(xi_in)=0 and f1(xi_out) = 1
	// When h_in=0, f1==f_F
	//
	vecd F(h.size());
	for (size_t i = 0; i < h.size(); ++i) {
		const double xi_LS2000 = h[i] / h.back();
		F[i] = F0 * oprel.f_F(xi_LS2000) * (1. - h.front() / h[i]) / (1. - h.front() / h.back());
	}
	return F;
}

double DiskStructureArguments::InitialFQuasistat::Coeff(double h_in, double h_out, const OpacityRelated &oprel) {
	odeint::runge_kutta_cash_karp54<double> stepper;
	const double x0 = h_in / (h_out - h_in);
	const double x1 = h_in / h_out;
	double integral = 0.;
	integrate_adaptive(
			stepper,
			[x0, x1, &oprel](const double &y, double &dydx, double x) {
				dydx = pow(oprel.f_F(x * (1. - x1) + x1) * x / (x * (1. - x1) + x1),
						   1. - oprel.m) * pow(x + x0, oprel.n);
			},
			integral, 0., 1., 0.01
	);
	return std::pow(h_out - h_in, oprel.n + 1.) * integral / ((1. - oprel.m) * oprel.D);
}

vecd DiskStructureArguments::InitialFGaussF::operator()(const vecd& h) const {
	vecd F(h.size());
	for (size_t i = 0; i < h.size(); ++i) {
		const double xi = (h[i] - h.front()) / (h.back() - h.front());
		F[i] = F0 * std::exp(-m::pow<2>(xi - gaussmu) / (2. * m::pow<2>(gausssigma)));
	}
	// Add some tiny accretion rate to avoid numerical problems
	if ((F[1] == 0.0) || (F[h.size()-1] == 0.0)) {
		for (size_t i = 0; i < h.size(); ++i) {
			const double xi = (h[i] - h.front()) / (h.back() - h.front());
			F[i] += 1e-100 * F0 * xi;
		}
	}
	return F;
}


constexpr const double SelfIrradiationArguments::default_Cirr;
constexpr const double SelfIrradiationArguments::default_irrindex;
constexpr const double SelfIrradiationArguments::default_Cirr_cold;
constexpr const double SelfIrradiationArguments::default_irrindex_cold;
constexpr const double SelfIrradiationArguments::default_height_to_radius_cold;
constexpr const char SelfIrradiationArguments::default_angular_dist_disk[];


constexpr const double FluxArguments::default_colourfactor;
constexpr const double FluxArguments::default_emin;
constexpr const double FluxArguments::default_emax;
constexpr const double FluxArguments::default_star_albedo;
constexpr const double FluxArguments::default_inclination;
constexpr const double FluxArguments::default_ephemeris_t0;


constexpr const double CalculationArguments::default_init_time;
constexpr const unsigned int CalculationArguments::default_Nx;
constexpr const unsigned int CalculationArguments::default_Nt_for_tau;
constexpr const char CalculationArguments::default_gridscale[];
constexpr const unsigned short CalculationArguments::default_starlod;
constexpr const int CalculationArguments::default_verb_level;
