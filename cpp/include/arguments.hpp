#ifndef _ARGUMENTS_HPP
#define _ARGUMENTS_HPP

#include <cmath>
#include <memory>
#include <optional>
#include <string>

#include <gsl_const_cgsm.h>
#include <constants.hpp>
#include <opacity_related.hpp>
#include <orbit.hpp>
#include <passband.hpp>
#include <util.hpp>
#include <unit_transformation.hpp>


class GeneralArguments {
public:
	constexpr static const char default_prefix[] = "freddi";
	constexpr static const char default_dir[] = ".";
	constexpr static const unsigned short default_output_precision = 12;
	constexpr static const unsigned int default_temp_sparsity_output = 1;
public:
	std::string prefix;
	std::string dir;
	unsigned short output_precision;
	unsigned int temp_sparsity_output;
	bool fulldata;
	bool stdout;
public:
	GeneralArguments(const std::string& prefix, const std::string& dir,
				  unsigned short output_precision,
				  unsigned int temp_sparsity_output,
				  bool fulldata,
				  bool stdout):
			prefix(prefix),
			dir(dir),
			output_precision(output_precision),
			temp_sparsity_output(temp_sparsity_output),
			fulldata(fulldata),
			stdout(stdout) {}
};


class BasicDiskBinaryArguments: public BlackHoleFunctions, public BinaryFunctions {
public:
	constexpr static const double default_alpha_to_alphacold = 10.0;
	constexpr static const double default_kerr = 0.0;
	constexpr static const double default_roche_lobe_fill = 1.0;
	constexpr static const double default_Topt = 0.0;
public:
	double alpha;
	double alphacold;
	double Mx;
	double kerr;
	double period;
	double Mopt;
	double roche_lobe_fill;
	double Topt;
	double rin;
	double rout;
	double risco;
protected:
	static inline double rinFromMxKerr(double Mx, double kerr) { return r_kerrISCO(Mx, kerr); }
	static inline double routFromMxMoptPeriod(double Mx, double Mopt, double period) {
		// 0.9 is approximation for r_max value from Paczyncki, 1977. He used grain model of accretting matter,
		// so his disk should be smaller than gas disk with pressure. So, probably, r_max is better approximation
		// than r_1 or r_2 for axial-symmetric gas disk. Also this factor is in agreement with Gilfanov & Arefiev, 2005
		return 0.9 * rocheLobeVolumeRadius(Mx, Mopt, period);
	}
	static inline double riscoFromMxKerr(double Mx, double kerr) { return r_kerrISCO(Mx, kerr); }
public:
	BasicDiskBinaryArguments(
			double alpha, std::optional<double> alphacold,
			double Mx, double kerr,
			double period,
			double Mopt, double roche_lobe_fill, double Topt,
			std::optional<double> rin, std::optional<double> rout, std::optional<double> risco
	):
			alpha(alpha), alphacold(alphacold ? *alphacold : alpha / default_alpha_to_alphacold),
			Mx(Mx), kerr(kerr),
			period(period),
			Mopt(Mopt),
			roche_lobe_fill(roche_lobe_fill),
			Topt(Topt),
			rin(rin ? *rin : rinFromMxKerr(Mx, kerr)),
			rout(rout ? *rout : routFromMxMoptPeriod(Mx, Mopt, period)),
			risco(risco ? *risco : riscoFromMxKerr(Mx, kerr)) {}
	BasicDiskBinaryArguments(const BasicDiskBinaryArguments&) = default;
	BasicDiskBinaryArguments(BasicDiskBinaryArguments&&) = default;
	inline double h(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx * r); }
	inline double omega(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx / r); }
};


class DiskStructureArguments {
protected:
	class InitialFFunction {
	protected:
		double F0;
		double Mdot0;
		double Mdisk0;
	public:
		InitialFFunction(double F0, double Mdot0, double Mdisk0):
				F0(F0), Mdot0(Mdot0), Mdisk0(Mdisk0) {}
		virtual ~InitialFFunction() = 0;
		virtual vecd operator()(const vecd& h) const = 0;
	};

	class InitialFPowerF: public InitialFFunction {
	protected:
		double powerorder;
	public:
		InitialFPowerF(double F0, double Mdot0, double Mdisk0, double powerorder):
				InitialFFunction(F0, Mdot0, Mdisk0),
				powerorder(powerorder) {}
		~InitialFPowerF() override = default;
		vecd operator()(const vecd& h) const override;
	};

	class InitialFLinearF: public InitialFPowerF {
	public:
		InitialFLinearF(double F0, double Mdot0, double Mdisk0):
				InitialFPowerF(F0, Mdot0, Mdisk0, 1.0) {}
	};

	class InitialFPowerSigma: public InitialFFunction {
	protected:
		double powerorder;
		std::shared_ptr<const OpacityRelated> oprel;
	public:
		InitialFPowerSigma(double F0, double Mdot0, double Mdisk0, double powerorder,
				std::shared_ptr<const OpacityRelated> oprel):
				InitialFFunction(F0, Mdot0, Mdisk0),
				powerorder(powerorder),
				oprel(oprel) {}
		~InitialFPowerSigma() override = default;
		vecd operator()(const vecd& h) const override;
	};

	class InitialFSineF: public InitialFFunction {
	public:
		using InitialFFunction::InitialFFunction;
		~InitialFSineF() override = default;
		vecd operator()(const vecd& h) const override;
	};

	class InitialFQuasistat: public InitialFFunction {
	protected:
		std::shared_ptr<const OpacityRelated> oprel;
	public:
		InitialFQuasistat(double F0, double Mdot0, double Mdisk0, const std::shared_ptr<const OpacityRelated>& oprel):
				InitialFFunction(F0, Mdot0, Mdisk0),
				oprel(oprel) {}
		~InitialFQuasistat() override = default;
		vecd operator()(const vecd& h) const override;
	};

	class InitialFGaussF: public InitialFFunction {
	protected:
		double gaussmu;
		double gausssigma;
	public:
		InitialFGaussF(double F0, double Mdot0, double Mdisk0, double gaussmu, double gausssigma):
				InitialFFunction(F0, Mdot0, Mdisk0),
				gaussmu(gaussmu), gausssigma(gausssigma) {}
		~InitialFGaussF() override = default;
		vecd operator()(const vecd& h) const override;
	};
public:
	constexpr static const char default_opacity[] = "Kramers";
	constexpr static const double default_Mdotout = 0.;
	constexpr static const char default_boundcond[] = "Teff";
	constexpr static const double default_Thot = 0.;
	constexpr static const double default_Tirr2Tvishot = 0.;
	constexpr static const char default_initialcond[] = "powerF";
	constexpr static const char default_wind[] = "no";
public:
	constexpr static const double mu = 0.62;
public:
	std::string opacity;
	std::shared_ptr<const OpacityRelated> oprel;
	double Mdotout;
	std::string boundcond;
	double Thot;
	double Tirr2Tvishot;
	double F0;
	double Mdisk0;
	double Mdot0;
	std::string initialcond;
	std::string wind;
	pard windparams;
protected:
	std::shared_ptr<InitialFFunction> initial_F_function;
public:
	DiskStructureArguments(
			const BasicDiskBinaryArguments &bdb_args,
			const std::string& opacity,
			double Mdotout,
			const std::string& boundcond, double Thot, double Tirr2Tvishot,
			const std::string& initialcond,
			std::optional<double> F0,
			std::optional<double> Mdisk0, std::optional<double> Mdot0,
			std::optional<double> powerorder,
			std::optional<double> gaussmu, std::optional<double> gausssigma,
			const std::string& wind, const pard& windparams);
	inline vecd initial_F(const vecd& h) const { return (*initial_F_function)(h); }
};


class SelfIrradiationArguments {
public:
	constexpr static const double default_Cirr = 0.0;
	constexpr static const double default_irrindex = 0.0;
	constexpr static const double default_Cirr_cold = 0.0;
	constexpr static const double default_irrindex_cold = 0.0;
	constexpr static const double default_height_to_radius_cold = 0.0;
	constexpr static const char default_angular_dist_disk[] = "plane";
public:
	double Cirr;
	double irrindex;
	double Cirr_cold;
	double irrindex_cold;
	double height_to_radius_cold;
	std::string angular_dist_disk;
public:
	SelfIrradiationArguments(
			double Cirr, double irrindex,
			double Cirr_cold, double irrindex_cold, double height_to_radius_cold,
			const std::string& angular_dist_disk):
			Cirr(Cirr), irrindex(irrindex),
			Cirr_cold(Cirr_cold), irrindex_cold(irrindex_cold), height_to_radius_cold(height_to_radius_cold),
			angular_dist_disk(angular_dist_disk) {}
};


class FluxArguments {
public:
	constexpr static const double default_colourfactor = 1.7;
	constexpr static const double default_emin = kevToHertz(1.);
	constexpr static const double default_emax = kevToHertz(12.);
	constexpr static const double default_star_albedo = 0.0;
	constexpr static const double default_inclination = 0.;  // degrees
	constexpr static const double default_ephemeris_t0 = 0.0;
public:
	double colourfactor;
	double emin;
	double emax;
	double star_albedo;
	double inclination;  // degrees
	double ephemeris_t0;
	double distance;
	bool cold_disk;
	bool star;
	vecd lambdas;
	std::vector<Passband> passbands;
public:
	FluxArguments(
			double colourfactor,
			double emin, double emax,
			double star_albedo,
			double inclination, double ephemeris_t0, double distance,
	        bool cold_disk, bool star,
	        const vecd& lambdas,
	        const std::vector<Passband>& passbands):
			colourfactor(colourfactor),
			emin(emin), emax(emax),
			star_albedo(star_albedo),
			inclination(inclination), ephemeris_t0(ephemeris_t0), distance(distance),
			cold_disk(cold_disk), star(star),
			lambdas(lambdas),
			passbands(passbands) {}
};


class CalculationArguments {
public:
	constexpr static const double default_init_time = 0;
	constexpr static const unsigned int default_Nx = 1000;
	constexpr static const unsigned int default_Nt_for_tau = 200;
	constexpr static const char default_gridscale[] = "log";
	constexpr static const unsigned short default_starlod = 3;
public:
	double init_time;
	double time;
	double tau;
	unsigned int Nx;
	std::string gridscale;
	unsigned short starlod = 3;
	double eps;
public:
	CalculationArguments(
			double inittime,
			double time, std::optional<double> tau,
			unsigned int Nx, const std::string& gridscale, const unsigned short starlod,
			double eps=1e-6):
			init_time(inittime),
			time(time),
			tau(tau ? *tau : time / default_Nt_for_tau),
			Nx(Nx), gridscale(gridscale), starlod(starlod),
			eps(eps) {}
};


class FreddiArguments {
public:
	std::shared_ptr<GeneralArguments> general;
	std::shared_ptr<BasicDiskBinaryArguments> basic;
	std::shared_ptr<DiskStructureArguments> disk;
	std::shared_ptr<SelfIrradiationArguments> irr;
	std::shared_ptr<FluxArguments> flux;
	std::shared_ptr<CalculationArguments> calc;
public:
	FreddiArguments() = default;
	FreddiArguments(
			GeneralArguments* general,
			BasicDiskBinaryArguments* basic,
			DiskStructureArguments* disk,
			SelfIrradiationArguments* irr,
			FluxArguments* flux,
			CalculationArguments* calc):
			general(general), basic(basic), disk(disk), irr(irr), flux(flux), calc(calc) {}
};

#endif // _ARGUMENTS_HPP
