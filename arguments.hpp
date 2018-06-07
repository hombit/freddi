#ifndef _ARGUMENTS_HPP
#define _ARGUMENTS_HPP

#include <cmath>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "gsl_const_cgsm.h"
#include "constants.hpp"
#include "orbit.hpp"


namespace po = boost::program_options;


class GeneralArguments {
public:
	constexpr static const char default_prefix[] = "freddi";
	constexpr static const char default_dir[] = "freddi";
public:
	std::string prefix;
	std::string dir;
	bool fulldata;
public:
	GeneralArguments(const char prefix[], const char dir[], bool fulldata=false) = default;
	static po::options_description description();
};


class MassTransformation {
public:
	static inline double gramToSun(const double mass_gram) { return mass_gram / GSL_CONST_CGSM_SOLAR_MASS; }
	static inline double sunToGram(const double mass_sun) { return mass_sun * GSL_CONST_CGSM_SOLAR_MASS; }
};


class LengthTransformation {
public:
	static inline double cmToRg(const double length_cm, const double mass_gram) {
		return length_cm * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
	}
	static inline double rgToCm(const double length_rg, const double mass_gram) {
		return length_rg * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / (GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * mass_gram);
	}
	static inline double kpcToCm(const double length_kpc) {
		return length_kpc * 1000. * GSL_CONST_CGSM_PARSEC;
	}
	static inline double cmToKpc(const double length_cm) {
		return length_cm / (1000. * GSL_CONST_CGSM_PARSEC);
	}
};


class TimeTransformation {
public:
	static inline double sToDay(double time_s) { return time_s / 86400.; }
	static inline double dayToS(double time_day) { return time_day * 86400.; }
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
	double alpha;
	double Mx;
	double kerr;
	double Mopt;
	double period;
	double rin;
	double rout;
public:
	BasicDiskBinaryArguments(double alpha,
							 double Mx, double kerr,
							 double Mopt, double period,
							 double rin, double rout) = default;
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
	std::string opacity;
	std::string boundcond;
	double Thot;
	std::string initialcond;
	double F0;
	double powerorder;
	double gaussmu;
	double gausssigma;
	double Mdisk0;
	double Mdot0;
public:
	DiskStructureArguments(const std::string& opacity,
						   const std::string& boundcond, double Thot,
						   const std::string& initialcond, double F0,
						   double powerorder, double gaussmu, double gausssigma,
						   double Mdisk0=-1, double Mdot0=-1) = default;
	static po::options_description description();
};


class SelfIrradiationArguments {
public:
	constexpr static const double default_Cirr = 0.;
	constexpr static const char default_irrfactortype[] = "const";
public:
	double Cirr;
	std::string irrfactortype;
	static po::options_description description();
};


class PhotonTransformation {
public:
	static inline double kevToHertz(const double energy_keV) {
		return energy_keV * 1000. * GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
	}
	static inline double hertzToKev(const double nu_hertz) {
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
	double colourfactor;
	double emin;
	double emax;
	double inclination;  // degrees
	double distance;
	std::vector<double> lambda;
public:
	static po::options_description description();
};


class CalculationArguments: public TimeTransformation {
public:
	constexpr static const double default_time = dayToS(50.);
	constexpr static const double default_tau = dayToS(0.25);
	constexpr static const unsigned int default_Nx = 1000;
	constexpr static const char default_gridscale[] = "log";
public:
	double time;
	double tau;
	unsigned int Nx;
	std::string gridscale;
	double eps;
public:
	CalculationArguments(double time, double tau, unsigned int Nx, const std::string& gridscale,
						 double eps=1e-6) = default;
	static po::options_description description();
};


class FreddiArguments {
public:
	GeneralArguments general;
	BasicDiskBinaryArguments basic;
	DiskStructureArguments disk;
	SelfIrradiationArguments irr;
	FluxArguments flux;
	CalculationArguments calc;
public:
	FreddiArguments(int argc, const char *argv[]);
};


#endif // _ARGUMENTS_HPP
