#ifndef FREDDI_NS_ARGUMENTS_HPP
#define FREDDI_NS_ARGUMENTS_HPP

#include <string>

#include <arguments.hpp>
#include <util.hpp>

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


#endif //FREDDI_NS_ARGUMENTS_HPP
