#ifndef _ARGUMENTS_HPP
#define _ARGUMENTS_HPP

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "gsl_const_cgsm.h"
#include "constants.hpp"
#include "opacity_related.hpp"
#include "orbit.hpp"


namespace po = boost::program_options;


class GeneralArguments {
public:
	constexpr static const char default_prefix[] = "freddi";
	constexpr static const char default_dir[] = ".";
public:
	const std::string prefix;
	const std::string dir;
	const bool fulldata;
public:
	GeneralArguments(const std::string& prefix, const std::string& dir, bool fulldata=false):
			prefix(prefix),
			dir(dir),
			fulldata(fulldata) {}
	GeneralArguments(const po::variables_map& vm);
	static po::options_description description();
};


class MassTransformation {
public:
	constexpr static inline double gramToSun(const double mass_gram) { return mass_gram / GSL_CONST_CGSM_SOLAR_MASS; }
	constexpr static inline double sunToGram(const double mass_sun) { return mass_sun * GSL_CONST_CGSM_SOLAR_MASS; }
};


class LengthTransformation {
public:
	constexpr static inline double cmToRg(const double length_cm, const double mass_gram) {
		return length_cm * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram);
	}
	constexpr static inline double rgToCm(const double length_rg, const double mass_gram) {
		return length_rg * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
	}
	constexpr static inline double kpcToCm(const double length_kpc) {
		return length_kpc * 1000. * GSL_CONST_CGSM_PARSEC;
	}
	constexpr static inline double cmToKpc(const double length_cm) {
		return length_cm / (1000. * GSL_CONST_CGSM_PARSEC);
	}
	static inline double sunToCm(const double length_solar_radius) {
		return length_solar_radius / solar_radius;
	}
	static inline double cmToSun(const double length_cm) {
		return length_cm * solar_radius;
	}
	constexpr static inline double angstromToCm(const double length_angstrom) {
		return length_angstrom * 1e-8;
	}
	constexpr static inline double cmToAngstrom(const double length_cm) {
		return length_cm * 1e8;
	}
};


class TimeTransformation {
public:
	constexpr static inline double sToDay(double time_s) { return time_s / 86400.; }
	constexpr static inline double dayToS(double time_day) { return time_day * 86400.; }
};


class BlackHoleFunctions: public LengthTransformation {
public:
	static double rISCORg(double kerr);
	static inline double rISCO(const double Mx, const double kerr) { return rgToCm(rISCORg(kerr), Mx); }
	static inline double accretionEfficiency(const double kerr) { return 1. - std::sqrt(1. - 2. / 3. / rISCORg(kerr)); }
};


class BinaryFunctions {
public:
	static double rocheLobeVolumeRadiusSemiaxis(double MxToMopt);
	static inline double semiaxis(const double total_mass, const double period) {
		return std::cbrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * total_mass * period*period / (4. * M_PI*M_PI));
	}
	static inline double semiaxis(const double Mx, const double Mopt, const double period) {
		return semiaxis(Mx + Mopt, period);
	}
	static inline double rocheLobeVolumeRadius(const double Mx, const double Mopt, const double period) {
		return rocheLobeVolumeRadiusSemiaxis(Mx / Mopt) * semiaxis(Mx, Mopt, period);
	}
};


class BasicDiskBinaryArguments:
		public MassTransformation, public BlackHoleFunctions, public BinaryFunctions, public TimeTransformation {
public:
	constexpr static const double default_alpha = 0.25;
	constexpr static const double default_Mx = sunToGram(5.);
	constexpr static const double default_kerr = 0.;
	constexpr static const double default_Mopt = sunToGram(0.5);
	constexpr static const double default_period = dayToS(0.25);
public:
	const double alpha;
	const double Mx;
	const double kerr;
	const double Mopt;
	const double period;
	const double rin;
	const double rout;
protected:
	static double rinInitializer(const po::variables_map& vm, double Mx, double kerr);
	static double routInitializer(const po::variables_map& vm, double Mx, double Mopt, double period);
public:
	BasicDiskBinaryArguments(
			double alpha,
			double Mx, double kerr,
			double Mopt, double period,
			double rin, double rout
	):
			alpha(alpha),
			Mx(Mx), kerr(kerr),
			Mopt(Mopt), period(period),
			rin(rin), rout(rout) {}
	static inline double rinFromMxKerr(double Mx, double kerr) { return rISCO(Mx, kerr); }
	static inline double routFromMxMoptPeriod(double Mx, double Mopt, double period) {
		return 0.8 * rocheLobeVolumeRadius(Mx, Mopt, period);
	}
	static inline BasicDiskBinaryArguments constructWithoutRinRout(const double alpha,
																   const double Mx, const double kerr,
																   const double Mopt, const double period) {
		return {alpha, Mx, kerr, Mopt, period,
				rinFromMxKerr(Mx, kerr), routFromMxMoptPeriod(Mx, Mopt, period)};
	}
	static inline BasicDiskBinaryArguments constructWithoutRin(const double alpha,
															   const double Mx, const double kerr,
															   const double Mopt, const double period,
															   const double rout) {
		return {alpha, Mx, kerr, Mopt, period,
				rinFromMxKerr(Mx, kerr), rout};
	}
	static inline BasicDiskBinaryArguments constructWithoutRout(const double alpha,
															    const double Mx, const double kerr,
															    const double Mopt, const double period,
															    const double rin) {
		return {alpha, Mx, kerr, Mopt, period,
				rin, routFromMxMoptPeriod(Mx, Mopt, period)};
	}
	BasicDiskBinaryArguments(const po::variables_map& vm);
	inline double h(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx * r); }
	inline double omega(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx / r); }
	static po::options_description description();
};


class DiskStructureArguments {
public:
	constexpr static const char default_opacity[] = "Kramers";
	constexpr static const char default_boundcond[] = "Teff";
	constexpr static const double default_Thot = 0.;
	constexpr static const char default_initialcond[] = "powerF";
	constexpr static const double default_F0 = 2e38;
	constexpr static const double default_powerorder = 6.;
	constexpr static const double default_gaussmu = 1.;
	constexpr static const double default_gausssigma = 0.25;
public:
	constexpr static const double mu = 0.62;
public:
	const std::string opacity;
	const std::string boundcond;
	const double Thot;
	const std::string initialcond;
	const double F0;
	const double powerorder;
	const double gaussmu;
	const double gausssigma;
	const double Mdisk0;
	const double Mdot0;
public:
	std::unique_ptr<const OpacityRelated> oprel;
protected:
	static double Mdisk0Initializer(const po::variables_map& vm);
	static double Mdot0Initializer(const po::variables_map& vm);
	static double F0Initializer(const po::variables_map& vm, const BasicDiskBinaryArguments& bdb_args);
public:
	DiskStructureArguments(
			const std::string& opacity,
			const std::string& boundcond, double Thot,
			const std::string& initialcond, double F0,
			double powerorder, double gaussmu, double gausssigma,
			double Mdisk0=-1, double Mdot0=-1):
			opacity(opacity),
			boundcond(boundcond), Thot(Thot),
			initialcond(initialcond), F0(F0),
			powerorder(powerorder), gaussmu(gaussmu), gausssigma(gausssigma),
			Mdisk0(Mdisk0), Mdot0(Mdot0) {}
	DiskStructureArguments(const po::variables_map& vm, const BasicDiskBinaryArguments& bdb_args);
	static po::options_description description();
};


class SelfIrradiationArguments {
public:
	constexpr static const double default_Cirr = 0.;
	constexpr static const char default_irrfactortype[] = "const";
public:
	const double Cirr;
	const std::string irrfactortype;
public:
	SelfIrradiationArguments(double Cirr, const std::string& irrfactortype):
			Cirr(Cirr), irrfactortype(irrfactortype) {}
	SelfIrradiationArguments(const po::variables_map& vm, const DiskStructureArguments& dsa_args);
	static po::options_description description();
};


class PhotonTransformation {
public:
	constexpr static inline double kevToHertz(const double energy_keV) {
		return energy_keV * 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
	}
	constexpr static inline double hertzToKev(const double nu_hertz) {
		return nu_hertz * GSL_CONST_CGSM_PLANCKS_CONSTANT_H / (1000. * GSL_CONST_CGSM_ELECTRON_VOLT);
	}
};


class FluxArguments: public PhotonTransformation, public LengthTransformation {
public:
	constexpr static const double default_colourfactor = 1.7;
	constexpr static const double default_emin = kevToHertz(1.);
	constexpr static const double default_emax = kevToHertz(12.);
	constexpr static const double default_inclination = 0.;  // degrees
	constexpr static const double default_distance = kpcToCm(10.);
public:
	const double colourfactor;
	const double emin;
	const double emax;
	const double inclination;  // degrees
	const double distance;
	const std::vector<double> lambdas;
protected:
	std::vector<double> lambdasInitializer(const po::variables_map& vm) const;
public:
	FluxArguments(
			double colourfactor,
			double emin, double emax,
			double inclination, double distance,
	        const std::vector<double>& lambdas):
			colourfactor(colourfactor),
			emin(emin), emax(emax),
			inclination(inclination), distance(distance),
			lambdas(lambdas) {}
	FluxArguments(const po::variables_map& vm);
	static po::options_description description();
};


class CalculationArguments: public TimeTransformation {
public:
	constexpr static const double default_time = dayToS(50.);
	constexpr static const double default_tau = dayToS(0.25);
	constexpr static const unsigned int default_Nx = 1000;
	constexpr static const char default_gridscale[] = "log";
public:
	const double time;
	const double tau;
	const unsigned int Nx;
	const std::string gridscale;
	const double eps;
public:
	CalculationArguments(
			double time, double tau, unsigned int Nx, const std::string& gridscale,
			double eps=1e-6):
			time(time), tau(tau), Nx(Nx), gridscale(gridscale), eps(eps) {}
	CalculationArguments(const po::variables_map& vm);
	static po::options_description description();
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
	FreddiArguments(
			GeneralArguments* general,
			BasicDiskBinaryArguments* basic,
			DiskStructureArguments* disk,
			SelfIrradiationArguments* irr,
			FluxArguments* flux,
			CalculationArguments* calc):
			general(general), basic(basic), disk(disk), irr(irr), flux(flux), calc(calc) {}
	FreddiArguments(const po::variables_map& vm);
	static po::options_description description();
};


po::variables_map parseArguments(int ac, char* av[]);


#endif // _ARGUMENTS_HPP
