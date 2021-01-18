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
		double alpha, const object& alphacold,
		double Mx, double kerr,
		double period,
		double Mopt, double roche_lobe_fill, double Topt,
		const object& rin, const object& rout, const object& risco);

boost::shared_ptr<DiskStructureArguments> make_disk_structure_arguments(
		const BasicDiskBinaryArguments& basic_disk_binary_arguments,
		const std::string& opacity,
		double Mdotout,
		const std::string& boundcond, double Thot, double Tirr2Tvishot,
		const std::string& initialcond,
		const object& F0, const object& Mdisk0, const object& Mdot0,
		const object& powerorder, const object& gaussmu, const object& gausssigma,
		const std::string& wind, const object& windparams);

boost::shared_ptr<SelfIrradiationArguments> make_self_irradiation_arguments(
		double Cirr, double irrindex,
		double Cirr_cold, double irrindex_cold, double height_to_radius_cold,
		const std::string& angular_dist_disk);

boost::shared_ptr<FluxArguments> make_flux_arguments(
			double colourfactor,
			double emin, double emax,
			double star_albedo,
			double inclination, double ephemeris_t0, double distance);

boost::shared_ptr<CalculationArguments> make_calculation_arguments(
		double inittime,
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
		const std::string& nsprop,
		const object& freqx, const object& Rx, double Bx, double hotspotarea,
		double epsilonAlfven, double inversebeta, double Rdead,
		const std::string& fptype, const object& fpparams_,
		const std::string& kappat_type, const object& kappat_params_,
		const std::string& ns_grav_redshift);

boost::shared_ptr<NeutronStarBasicDiskBinaryArguments> make_neutron_star_basic_disk_binary_arguments(
		const NeutronStarArguments& ns_args,
		double alpha, const object& alphacold,
		double Mx, double kerr,
		double period,
		double Mopt, double roche_lobe_fill, double Topt,
		const object& rin, const object& rout, const object& risco);

boost::shared_ptr<NeutronStarSelfIrradiationArguments> make_neutron_star_self_irradiation_arguments(
		double Cirr, double irrindex,
		double Cirr_cold, double irrindex_cold, double default_height_to_radius_cold,
		const std::string& angular_dist_disk, const std::string& angular_dist_ns);

boost::shared_ptr<FreddiNeutronStarArguments> make_freddi_neutron_star_arguments(
		const GeneralArguments& general,
		const NeutronStarBasicDiskBinaryArguments& basic,
		const DiskStructureArguments& disk,
		const NeutronStarSelfIrradiationArguments& irr,
		const FluxArguments& flux,
		const CalculationArguments& calc,
		const NeutronStarArguments& ns);


#endif //FREDDI_PYTHON_ARGUMENTS_HPP
