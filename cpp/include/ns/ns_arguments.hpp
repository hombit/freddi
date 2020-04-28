#ifndef FREDDI_NS_ARGUMENTS_HPP
#define FREDDI_NS_ARGUMENTS_HPP

#include <optional>
#include <map>
#include <string>

#include <arguments.hpp>
#include <util.hpp>

namespace SibgatullinSunyaev2000Geometry {
	static double radiusNS(double freqx);
	static double radiusISCO(double reqx);
};

class NeutronStarArguments {
public:
	constexpr static const char default_nsprop[] = "dummy";
	constexpr static const double default_hotspotarea = 1.;
	constexpr static const double default_epsilonAlfven = 1.;
	constexpr static const double default_inversebeta = 0.;
	constexpr static const double default_Rdead = 0.;
	constexpr static const char default_fptype[] = "no-outflow";
	constexpr static const char default_kappat_type[] = "const";
	constexpr static const double default_kappat_value = 1.0 / 3.0;
	constexpr static const char default_ns_grav_redshift[] = "off";
protected:
	constexpr static const double default_Rx_dummy = 1e6;
	constexpr static const double default_freqx_dummy = 0.;
public:
	std::string nsprop;
	double freqx;
	double Rx;
	double Bx;
	double hotspotarea;
	double epsilonAlfven;
	double inversebeta;
	double Rdead;
	std::string fptype;
	pard fpparams;
	std::string kappat_type;
	pard kappat_params;
	std::string ns_grav_redshift;
protected:
	static double initializeFreqx(const std::string& nsprop);
	static double initializeRx(const std::string& nsprop, std::optional<double> freqx);
public:
	NeutronStarArguments(
			const std::string& nsprop,
			std::optional<double> freqx,
			std::optional<double> Rx,
			double Bx, double hotspotarea,
			double epsilonAlfven, double inversebeta, double Rdead,
			const std::string& fptype, const pard& fpparams,
			const std::string& kappat_type, const pard& kappat_params,
			const std::string& ns_grav_redshift):
			nsprop(nsprop),
			freqx(freqx ? *freqx : initializeFreqx(nsprop)),
			Rx(Rx ? *Rx : initializeRx(nsprop, freqx)),
			Bx(Bx), hotspotarea(hotspotarea),
			epsilonAlfven(epsilonAlfven), inversebeta(inversebeta), Rdead(Rdead),
			fptype(fptype), fpparams(fpparams),
			kappat_type(kappat_type), kappat_params(kappat_params),
			ns_grav_redshift(ns_grav_redshift) {}
};


class NeutronStarBasicDiskBinaryArguments: public BasicDiskBinaryArguments {
protected:
	static std::optional<double> initializeRisco(const NeutronStarArguments& ns_args, std::optional<double> risco);
public:
	NeutronStarBasicDiskBinaryArguments(
			const NeutronStarArguments& ns_args,
			double alpha_, std::optional<double> alphacold_,
			double Mx_, double kerr_,
			double period_,
			double Mopt_, double roche_lobe_fill_, double Topt_,
			std::optional<double> rin_, std::optional<double> rout_, std::optional<double> risco_
	);
};


class NeutronStarSelfIrradiationArguments: public SelfIrradiationArguments {
public:
	constexpr static const char default_angular_dist_ns[] = "isotropic";
public:
	std::string angular_dist_ns;
public:
	NeutronStarSelfIrradiationArguments(
			double Cirr, double irrindex,
			double Cirr_cold, double irrindex_cold, double height_to_radius_cold,
			const std::string& angular_dist_disk, const std::string& angular_dist_ns):
			SelfIrradiationArguments(Cirr, irrindex, Cirr_cold,
					irrindex_cold, height_to_radius_cold,
					angular_dist_disk),
			angular_dist_ns(angular_dist_ns) {}
};


class FreddiNeutronStarArguments: public FreddiArguments {
public:
	std::shared_ptr<NeutronStarArguments> ns;
	std::shared_ptr<NeutronStarSelfIrradiationArguments> irr_ns;
public:
	FreddiNeutronStarArguments() = default;
	FreddiNeutronStarArguments(
			GeneralArguments* general_,
			NeutronStarBasicDiskBinaryArguments* basic_,
			DiskStructureArguments* disk_,
			NeutronStarSelfIrradiationArguments* irr_,
			FluxArguments* flux_,
			CalculationArguments* calc_,
			NeutronStarArguments* ns_):
			ns(ns_), irr_ns(irr_) {
		general.reset(general_);
		basic.reset(basic_);
		disk.reset(disk_);
		irr = irr_ns;
		flux.reset(flux_);
		calc.reset(calc_);
	}
};


#endif //FREDDI_NS_ARGUMENTS_HPP
