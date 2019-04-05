#ifndef FREDDI_PYTHON_ARGUMENTS_HPP
#define FREDDI_PYTHON_ARGUMENTS_HPP

#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>

#include <arguments.hpp>

using object = boost::python::object;

void wrap_arguments();

boost::shared_ptr<GeneralArguments> make_general_arguments();

boost::shared_ptr<BasicDiskBinaryArguments> make_basic_disk_binary_arguments(
		double alpha,
		double Mx, double kerr,
		double Mopt, double period,
		const object& rin, const object& rout);

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot,
		const std::string& initialcond, double F0,
		double powerorder, double gaussmu, double gausssigma,
		const object& Mdisk0_, const object& Mdot0_,
		const std::string& wind, const object& windparams_);

boost::shared_ptr<SelfIrradiationArguments> make_self_irradiation_arguments(double Cirr, const std::string& irrfactortype);

boost::shared_ptr<FluxArguments> make_flux_arguments(
			double colourfactor,
			double emin, double emax,
			double inclination, double distance,
			const object& lambdas_);

boost::shared_ptr<CalculationArguments> make_calculation_arguments(
		double time, double tau, unsigned int Nx, const std::string& gridscale,
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
		const std::string& fptype);

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const FreddiArguments& freddi_args, const NeutronStarArguments& ns);


#endif //FREDDI_PYTHON_ARGUMENTS_HPP
