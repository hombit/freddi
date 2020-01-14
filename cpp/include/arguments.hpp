#ifndef _ARGUMENTS_HPP
#define _ARGUMENTS_HPP

#include <cmath>
#include <memory>
#include <string>

#include "gsl_const_cgsm.h"
#include "constants.hpp"
#include "opacity_related.hpp"
#include "orbit.hpp"
#include "passband.hpp"
#include "util.hpp"
#include "unit_transformation.hpp"


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
};


class BlackHoleFunctions {
public:
	static inline double rISCO(const double Mx, const double kerr) { return rgToCm(rISCORg(kerr), Mx); }
};


class BinaryFunctions {
public:
	static double rocheLobeVolumeRadiusSemiaxis(double MxToMopt);
	static inline double semiaxis(const double total_mass, const double period) {
		return std::cbrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * total_mass * m::pow<2>(period) / (4. * m::pow<2>(M_PI)));
	}
	static inline double semiaxis(const double Mx, const double Mopt, const double period) {
		return semiaxis(Mx + Mopt, period);
	}
	static inline double rocheLobeVolumeRadius(const double Mx, const double Mopt, const double period) {
		return rocheLobeVolumeRadiusSemiaxis(Mx / Mopt) * semiaxis(Mx, Mopt, period);
	}
};


class BasicDiskBinaryArguments: public BlackHoleFunctions, public BinaryFunctions {
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
	BasicDiskBinaryArguments(const BasicDiskBinaryArguments&) = default;
	BasicDiskBinaryArguments(BasicDiskBinaryArguments&&) = default;
	static inline double rinFromMxKerr(double Mx, double kerr) { return rISCO(Mx, kerr); }
	static inline double routFromMxMoptPeriod(double Mx, double Mopt, double period) {
		// 0.9 is approximation for r_max value from Paczyncki, 1977. He used grain model of accretting matter,
		// so his disk should be smaller than gas disk with pressure. So, probably, r_max is better approximation
		// than r_1 or r_2 for axial-symmetric gas disk. Also this factor is in agreement with Gilfanov & Arefiev, 2005
		return 0.9 * rocheLobeVolumeRadius(Mx, Mopt, period);
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
	inline double h(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx * r); }
	inline double omega(const double r) const { return std::sqrt(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx / r); }
};


class DiskStructureArguments {
public:
	constexpr static const char default_opacity[] = "Kramers";
	constexpr static const double default_Mdotout = 0.;
	constexpr static const char default_boundcond[] = "Teff";
	constexpr static const double default_Thot = 0.;
	constexpr static const char default_initialcond[] = "powerF";
	constexpr static const double default_F0 = 2e38;
	constexpr static const double default_powerorder = 6.;
	constexpr static const double default_gaussmu = 1.;
	constexpr static const double default_gausssigma = 0.25;
	constexpr static const char default_wind[] = "no";
public:
	constexpr static const double mu = 0.62;
public:
	const std::string opacity;
	const double Mdotout;
	const std::string boundcond;
	const double Thot;
	const std::string initialcond;
	const double powerorder;
	const double gaussmu;
	const double gausssigma;
	const std::string wind;
	const pard windparams;
	const double Mdisk0;
	const double Mdot0;
	const bool is_Mdisk0_specified;
	const bool is_Mdot0_specified;
	std::shared_ptr<const OpacityRelated> oprel;
	const double F0;
protected:
	double F0Initializer(double F0_, const BasicDiskBinaryArguments& bdb_args);
public:
	DiskStructureArguments(
			const BasicDiskBinaryArguments &bdb_args,
			const std::string& opacity,
			double Mdotout,
			const std::string& boundcond, double Thot,
			const std::string& initialcond, double F0,
			double powerorder, double gaussmu, double gausssigma,
			bool is_Mdisk0_specified, bool is_Mdot0_specified,
			double Mdisk0, double Mdot0,
			const std::string& wind, const pard& windparams);
	DiskStructureArguments(const DiskStructureArguments&) = default;
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
};


class FluxArguments {
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
	const vecd lambdas;
	const std::vector<Passband> passbands;
public:
	FluxArguments(
			double colourfactor,
			double emin, double emax,
			double inclination, double distance,
	        const vecd& lambdas,
	        const std::vector<Passband> passbands):
			colourfactor(colourfactor),
			emin(emin), emax(emax),
			inclination(inclination), distance(distance),
			lambdas(lambdas),
			passbands(passbands) {}
};


class CalculationArguments {
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


class NeutronStarArguments {
public:
	constexpr static const double default_Rx = 1e6;
	constexpr static const double default_freqx = 0.;
	constexpr static const double default_Bx = 0.;
	constexpr static const double default_hotspotarea = 1.;
	constexpr static const double default_epsilonAlfven = 1.;
	constexpr static const double default_inversebeta = 0.;
	constexpr static const double default_Rdead = 0.;
	constexpr static const char default_fptype[] = "no-outflow";
	const pard default_fpparams = {};
public:
	const double Rx;
	const double freqx;
	const double Bx;
	const double hotspotarea;
	const double epsilonAlfven;
	const double inversebeta;
	const double Rdead;
	const std::string fptype;
	const pard fpparams;
public:
	NeutronStarArguments(
			double Rx, double freqx, double Bx, double hotspotarea,
			double epsilonAlfven, double inversebeta, double Rdead,
			const std::string& fptype, const pard& fpparams):
			Rx(Rx), freqx(freqx), Bx(Bx), hotspotarea(hotspotarea),
			epsilonAlfven(epsilonAlfven), inversebeta(inversebeta), Rdead(Rdead),
			fptype(fptype), fpparams(fpparams) {}
};


class FreddiNeutronStarArguments: public FreddiArguments {
public:
	std::shared_ptr<NeutronStarArguments> ns;
public:
	FreddiNeutronStarArguments() = default;
	FreddiNeutronStarArguments(const FreddiArguments& freddi_args, NeutronStarArguments* ns):
			FreddiArguments(freddi_args), ns(ns) {}
	FreddiNeutronStarArguments(
			GeneralArguments* general,
			BasicDiskBinaryArguments* basic,
			DiskStructureArguments* disk,
			SelfIrradiationArguments* irr,
			FluxArguments* flux,
			CalculationArguments* calc,
			NeutronStarArguments* ns):
			FreddiArguments(general, basic, disk, irr, flux, calc), ns(ns) {}
};


#endif // _ARGUMENTS_HPP
