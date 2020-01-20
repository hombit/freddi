#include <boost/numeric/odeint.hpp>

#include "arguments.hpp"
#include "util.hpp"


namespace odeint = boost::numeric::odeint;


constexpr const char GeneralArguments::default_prefix[];
constexpr const char GeneralArguments::default_dir[];


constexpr const double BasicDiskBinaryArguments::default_alpha;
constexpr const double BasicDiskBinaryArguments::default_Mx;
constexpr const double BasicDiskBinaryArguments::default_kerr;
constexpr const double BasicDiskBinaryArguments::default_Mopt;
constexpr const double BasicDiskBinaryArguments::default_Topt;
constexpr const double BasicDiskBinaryArguments::default_period;


double BinaryFunctions::rocheLobeVolumeRadiusSemiaxis(const double mass_ratio) { // Eggleton, P. P. 1983, ApJ, 268, 368
	const double q = std::cbrt(mass_ratio);
	return 0.49 * m::pow<2>(q) / (0.6 * m::pow<2>(q) + std::log(1. + q));
}

constexpr const char DiskStructureArguments::default_opacity[];
constexpr const double DiskStructureArguments::default_Mdotout;
constexpr const char DiskStructureArguments::default_initialcond[];
constexpr const double DiskStructureArguments::default_Thot;
constexpr const char DiskStructureArguments::default_boundcond[];
constexpr const double DiskStructureArguments::default_F0;
constexpr const double DiskStructureArguments::default_powerorder;
constexpr const double DiskStructureArguments::default_gaussmu;
constexpr const double DiskStructureArguments::default_gausssigma;
constexpr const double DiskStructureArguments::mu;
constexpr const char DiskStructureArguments::default_wind[];

DiskStructureArguments::DiskStructureArguments(
	const BasicDiskBinaryArguments &bdb_args,
	const std::string &opacity,
	const double Mdotout,
	const std::string &boundcond, const double Thot,
	const std::string &initialcond,
	const double F0,
	const double powerorder, const double gaussmu, const double gausssigma,
	const bool is_Mdisk0_specified, const bool is_Mdot0_specified,
	const double Mdisk0, const double Mdot0,
	const std::string& wind, const pard& windparams):
		opacity(opacity),
		oprel(new OpacityRelated(opacity, bdb_args.Mx, bdb_args.alpha, mu)),
		Mdotout(Mdotout),
		boundcond(boundcond),
		Thot(Thot),
		initialcond(initialcond),
		powerorder(powerorder),
		gaussmu(gaussmu),
		gausssigma(gausssigma),
		wind(wind),
		windparams(windparams),
		is_Mdisk0_specified(is_Mdisk0_specified),
		is_Mdot0_specified(is_Mdot0_specified),
		Mdisk0(Mdisk0),
		Mdot0(Mdot0),
		F0(F0Initializer(F0, bdb_args)) {}


double DiskStructureArguments::F0Initializer(double F0_, const BasicDiskBinaryArguments& bdb_args) {
	const double h_in = bdb_args.h(bdb_args.rin);
	const double h_out = bdb_args.h(bdb_args.rout);

	if (initialcond == "powerF" || initialcond == "power") {
		if (is_Mdot0_specified) {
			throw std::runtime_error("Set Mdisk0 or F0 instead of Mdot0 if initialcond=powerF !");
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = (1. - oprel->m) * powerorder;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a,b,x0]( const double &y, double &dydx, double x ){
						dydx = pow(x, a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0_ = pow(Mdisk0 * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral,
					 1. / (1. - oprel->m));
			return F0_;
		}
		return F0_;
	}
	if (initialcond == "powerSigma") {
		if (is_Mdot0_specified) {
			throw std::runtime_error("Set Mdisk0 or F0 instead of Mdot0 if initialcond=powerSigma !");
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = powerorder;
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
			F0_ = pow(Mdisk0 * (1. - oprel->m) * oprel->D * pow(h_out, 3. - oprel->n) / m::pow<4>(h_out - h_in) / integral,
					  1. / (1. - oprel->m));
			return F0_;
		}
		return F0_;
	}
	if (initialcond == "sinusF" || initialcond == "sinus") {
		if (is_Mdot0_specified) {
			F0_ = Mdot0 * (h_out - h_in) * 2./M_PI;
			return F0_;
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel->m;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a, b, x0](const double &y, double &dydx, double x) {
						dydx = pow(sin(x * M_PI_2), a) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0_ = pow(Mdisk0 * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral,
					  1. / (1. - oprel->m));
			return F0_;
		}
		return F0_;
	}
	if (initialcond == "gaussF") {
		if (gaussmu <= 0. or gaussmu > 1.){
			throw std::runtime_error("gaussmu value should be in [0..1] !");
		}
		if (is_Mdot0_specified) {
			throw std::runtime_error("Setting Mdot if initialcond=gaussF produces unstable results and it isn't motivated physically. Set F0 or Mdisk0 instead");
			//F0_ = Mdot_in * (h_out - h_in) * m::pow<2>(gauss_sigma) / gauss_mu * exp( m::pow<2>(gauss_mu) / (2. * m::pow<2>(gauss_sigma)) );
		}
		if (is_Mdisk0_specified){
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double a = 1. - oprel->m;
			const double b = oprel->n;
			const double x0 = h_in / (h_out - h_in);
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[a,b,x0,this]( const double &y, double &dydx, double x ){
						dydx = exp( -m::pow<2>(x - this->gaussmu) * a / (2. * m::pow<2>(this->gausssigma)) ) * pow(x + x0, b);
					},
					integral, 0., 1., 0.01
			);
			F0_ = pow(Mdisk0 * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral,
					  1. / (1. - oprel->m));
			return F0_;
		}
		return F0_;
	}
	if (initialcond == "quasistat") {
		if (is_Mdot0_specified) {
			F0_ = Mdot0 * (h_out - h_in) / h_out * h_in / oprel->f_F(h_in/h_out);
			return F0_;
		}
		if (is_Mdisk0_specified) {
			odeint::runge_kutta_cash_karp54<double> stepper;
			const double x0 = h_in / (h_out - h_in);
			const double x1 = h_in / h_out;
			double integral = 0.;
			integrate_adaptive(
					stepper,
					[x0,x1,this]( const double &y, double &dydx, double x ){
						dydx = pow(this->oprel->f_F(x * (1. - x1) + x1) * x / (x * (1. - x1) + x1), 1. - this->oprel->m) * pow(x + x0, this->oprel->n);
					},
					integral, 0., 1., 0.01
			);
			F0_ = pow(Mdisk0 * (1. - oprel->m) * oprel->D / pow(h_out - h_in, oprel->n + 1.) / integral, 1. / (1. - oprel->m));
			return F0_;
		}
		return F0_;
	}
	throw std::runtime_error("Invalid value of initialcond");
}


constexpr const double SelfIrradiationArguments::default_Cirr;
constexpr const char SelfIrradiationArguments::default_irrfactortype[];


constexpr const double FluxArguments::default_colourfactor;
constexpr const double FluxArguments::default_emin;
constexpr const double FluxArguments::default_emax;
constexpr const double FluxArguments::default_inclination;
constexpr const double FluxArguments::default_distance;
constexpr const bool FluxArguments::default_cold_disk;


constexpr const double CalculationArguments::default_time;
constexpr const double CalculationArguments::default_tau;
constexpr const unsigned int CalculationArguments::default_Nx;
constexpr const char CalculationArguments::default_gridscale[];
constexpr const unsigned short CalculationArguments::default_starlod;
