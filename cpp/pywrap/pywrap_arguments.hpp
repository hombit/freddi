#ifndef FREDDI_PYTHON_ARGUMENTS_HPP
#define FREDDI_PYTHON_ARGUMENTS_HPP

#include <optional>
#include <string>

#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>

#include <arguments.hpp>
#include <ns/ns_arguments.hpp>

using object = boost::python::object;

void wrap_arguments();

boost::shared_ptr<GeneralArguments> make_general_arguments();

boost::shared_ptr<BasicDiskBinaryArguments> make_basic_disk_binary_arguments(
		double alpha,
		double Mx, double kerr,
		double period,
		double Mopt, const object& Ropt, double Topt,
		const object& rin, const object& rout, const object& risco);

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot,
		const std::string& initialcond,
		const object& F0, const object& Mdisk0, const object& Mdot0,
		const object& powerorder, const object& gaussmu, const object& gausssigma,
		const std::string& wind, const object& windparams);

boost::shared_ptr<SelfIrradiationArguments> make_self_irradiation_arguments(double Cirr, const std::string& irrfactortype);

boost::shared_ptr<FluxArguments> make_flux_arguments(
			double colourfactor,
			double emin, double emax,
			double star_albedo,
			double inclination, double distance);

boost::shared_ptr<CalculationArguments> make_calculation_arguments(
		double time, const object& tau,
		unsigned int Nx, const std::string& gridscale, unsigned short starlod,
		const object& eps);

boost::shared_ptr<FreddiArguments> make_freddi_arguments(
		const GeneralArguments& general,
		const BasicDiskBinaryArguments& basic,
		const DiskStructureArguments& disk,
		const SelfIrradiationArguments& irr,
		const FluxArguments& flux,
		const CalculationArguments& calc);

boost::shared_ptr<NeutronStarArguments> make_neutron_star_arguments(
		double Rx, double freqx, double Bx, double hotspotarea,
		double epsilonAlfven, double inversebeta, double Rdead,
		const std::string& fptype, const object& fpparams_);

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const FreddiArguments& freddi_args, const NeutronStarArguments& ns);


#endif //FREDDI_PYTHON_ARGUMENTS_HPP
