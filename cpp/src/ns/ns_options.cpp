#include <ns/ns_options.hpp>

NeutronStarOptions::NeutronStarOptions(const po::variables_map &vm):
		NeutronStarArguments(
				vm["nsprop"].as<std::string>(),
				varToOpt<double>(vm, "freqx"),
				varToOpt<double>(vm, "Rx"),
				vm["Bx"].as<double>(),
				vm["hotspotarea"].as<double>(),
				vm["epsilonAlfven"].as<double>(),
				vm["inversebeta"].as<double>(),
				vm["Rdead"].as<double>(),
				vm["fptype"].as<std::string>(),
				multitoken_string_option_to_map(vm, "fpparams", ":", {}),
				vm["kappattype"].as<std::string>(),
				multitoken_string_option_to_map(vm, "kappatparams", ":", default_kappat_params.at(vm["kappattype"].as<std::string>())),
				vm["nsgravredshift"].as<std::string>()) {}


po::options_description NeutronStarOptions::description() {
	po::options_description od("Parameters of accreting neutron star");
	od.add_options()
			( "nsprop", po::value<std::string>()->default_value(default_nsprop), "Neutron star geometry and radiation properties name and specifies default values of --Rx, --Risco and --freqx\n\n"
																				 "Values:\n"
																				 "  dummy: NS radiation efficiency is R_g * (1 / R_x - 1 / 2R_in), default --freqx is 0, default Rx is 1e6, default Risco is Kerr value\n"
																				 "  newt: NS radiation efficiency is a functions of NS frequency, that's why --freqx must be specified explicitly\n"
																				 "  sibgatullin-sunyaev2000: NS radiation efficiency, and default values of Rx and Risco are functions of NS frequency, that's why --freqx must be specified explicitly")
			( "freqx", po::value<double>(), "Accretor rotation frequency, Hz. This parameter is not linked to --kerr, agree them yourself" )
			( "Rx", po::value<double>(), "Accretor radius, cm" )
			( "Bx", po::value<double>()->required(), "Accretor polar magnetic induction, G" )
			( "hotspotarea", po::value<double>()->default_value(default_hotspotarea), "??? Relative area of hot spot on the accretor" )
			( "epsilonAlfven", po::value<double>()->default_value(default_epsilonAlfven), "Factor in Alfven radius formula" )
			( "inversebeta", po::value<double>()->default_value(default_inversebeta), "???" )
			( "Rdead", po::value<double>()->default_value(default_Rdead), "Maximum inner radius of the disk that can be obtained, it characterises minimum torque in the dead disk, cm" )
			( "fptype", po::value<std::string>()->default_value(default_fptype), "Accretor Mdot fraction mode: no-outflow, propeller, corotation-block, eksi-kultu2010, romanova2018, geometrical" )
			( "fpparams", po::value<std::vector<std::string>>()->multitoken()->composing(), "Accretor Mdot fraction parameters, specific for each fptype. Format is name:value. Examples: 1) for geometrical chi:15; 2) for romanova2018 par1:0.15 par2:0.92" )
			( "kappattype", po::value<std::string>()->default_value(default_kappat_type), "kappa_t describes how strong is interaction between neutron star magnitosphere and disk, magnetic torque is kappa_t(R) * mu^2 / R^3. This parameter describes type of radial destribution of this parameter\n\n"
																						  "Values:\n"
																						  "  const: doesn't depend on radius, default value is 1/3, --kappatype name: 'value'\n"
																						  "  corstep: kappa_t takes one value inside corotation radius, and another outside, default values are 1/3, --kappatype names: 'in', 'out'\n"
																						  "  romanova2018: similar to corstep option, but the outside value is reduced by the portion taken away by wind (see Table 2 of Romanova+2018,NewA,62,94), --kappatype names: 'in', 'out'")
			( "kappatparams", po::value<std::vector<std::string>>()->multitoken()->composing(), "Parameters of kappa_t radial distribution, see --kappattype for details. Format is name:value" )
			( "nsgravredshift", po::value<std::string>()->default_value(default_ns_grav_redshift), "Neutron star gravitational redshift type.\n"
																						           "Values:\n"
								                                                                   "  off: gravitational redshift doesn't taken into account\n"
																		                           "  on: redshift is (1 - R_sch / Rx), where R_sch = 2GM/c^2" )
			;
	return od;
}


NeutronStarBasicDiskBinaryOptions::NeutronStarBasicDiskBinaryOptions(const po::variables_map &vm, const NeutronStarArguments& args_ns):
		NeutronStarBasicDiskBinaryArguments(
				args_ns,
				vm["alpha"].as<double>(),
				sunToGram(vm["Mx"].as<double>()),
				vm["kerr"].as<double>(),
				dayToS(vm["period"].as<double>()),
				sunToGram(vm["Mopt"].as<double>()),
				vm["rochelobefill"].as<double>(),
				vm["Topt"].as<double>(),
				BasicDiskBinaryOptions::rinInitializer(vm),
				BasicDiskBinaryOptions::routInitializer(vm),
				BasicDiskBinaryOptions::riscoInitializer(vm)) {}

po::options_description NeutronStarBasicDiskBinaryOptions::description() {
	return BasicDiskBinaryOptions::description();
}


NeutronStarSelfIrradiationOptions::NeutronStarSelfIrradiationOptions(const po::variables_map& vm, const DiskStructureArguments& dsa_args):
		NeutronStarSelfIrradiationArguments(
				vm["Cirr"].as<double>(),
				vm["irrindex"].as<double>(),
				vm["Cirrcold"].as<double>(),
				vm["irrindexcold"].as<double>(),
				vm["angulardistdisk"].as<std::string>(),
				vm["angulardistns"].as<std::string>()) {
	if (Cirr <= 0. && dsa_args.boundcond == "Tirr") {
		throw po::error("Set positive --Cirr when --boundcond=Tirr");
	}
}

po::options_description NeutronStarSelfIrradiationOptions::description() {
	auto od = SelfIrradiationOptions::description();
	od.add_options()
			( "angulardistns", po::value<std::string>()->default_value(default_angular_dist_ns), "Angular distribution type of the neutron star X-ray radiation. Values: isotropic, plane" )
			;
	return od;
}


FreddiNeutronStarOptions::FreddiNeutronStarOptions(const po::variables_map &vm) {
	ns.reset(new NeutronStarOptions(vm));

	general.reset(new GeneralOptions(vm));
	basic.reset(new NeutronStarBasicDiskBinaryOptions(vm, *ns));
	disk.reset(new DiskStructureOptions(vm, *basic));
	irr_ns.reset(new NeutronStarSelfIrradiationOptions(vm, *disk));
	irr = irr_ns;
	flux.reset(new FluxOptions(vm));
	calc.reset(new CalculationOptions(vm));
}

po::options_description FreddiNeutronStarOptions::description() {
	po::options_description desc("Freddi NS: numerical calculation of accretion disk evolution");
	desc.add(GeneralOptions::description());
	desc.add(NeutronStarBasicDiskBinaryOptions::description());
	desc.add(DiskStructureOptions::description());
	desc.add(NeutronStarOptions::description());
	desc.add(NeutronStarSelfIrradiationOptions::description());
	desc.add(FluxOptions::description());
	desc.add(CalculationOptions::description());
	return desc;
}

