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
				fpparamsInitializer(vm)) {}

pard NeutronStarOptions::fpparamsInitializer(const po::variables_map& vm) {
	pard m;
	if (vm.count("fpparams") == 0) {
		return m;
	}
	const std::vector<std::string > tokens(vm["fpparams"].as<std::vector<std::string> >());
	for (const auto& token : tokens) {
		std::vector<std::string> parts;
		boost::split(parts, token, boost::is_any_of(":"));
		if (parts.size() != 2) {
			throw po::validation_error(po::validation_error::invalid_option_value);
		}
		m[parts[0]] = std::stod(parts[1]);
	}
	return m;
}

po::options_description NeutronStarOptions::description() {
	po::options_description od("Parameters of accreting neutron star");
	od.add_options()
			( "nsprop", po::value<std::string>()->default_value(default_nsprop), "Neutron star geometry and radiation properties name and specifies default values of --Rx, --Risco and --freqx\n\n"
																				 "Values:\n"
																				 "  dummy: NS radiation efficiency is R_g * (1 / R_x - 1 / 2R_in), default --freqx is 0, default Rx is 1e6, default Risco is Kerr value\n"
																				 "  sibgatullin-sunyaev2000: NS radiation efficiency, and default values of Rx and Risco are functions of NS frequency, that's why --freqx must be specified explicitly")
			( "freqx", po::value<double>(), "Accretor rotation frequency, Hz. This parameter is not linked to --kerr, agree them yourself" )
			( "Rx", po::value<double>(), "Accretor radius, cm" )
			( "Bx", po::value<double>()->required(), "Accretor polar magnetic induction, G" )
			( "hotspotarea", po::value<double>()->default_value(default_hotspotarea), "??? Relative area of hot spot on the accretor" )
			( "epsilonAlfven", po::value<double>()->default_value(default_epsilonAlfven), "Factor in Alfven radius formula" )
			( "inversebeta", po::value<double>()->default_value(default_inversebeta), "???" )
			( "Rdead", po::value<double>()->default_value(default_Rdead), "Maximum inner radius of the disk that can be obtained, it characterises minimum torque in the dead disk, cm" )
			( "fptype", po::value<std::string>()->default_value(default_fptype), "Accretor Mdot fraction mode: no-outflow, propeller, corotation-block, eksi-kultu2010, romanova2018, geometrical" )
			( "fpparams", po::value<std::vector<std::string> >()->multitoken(), "Accretor Mdot fraction parameters, specific for each fptype. Format is name:value. Examples: 1) for geometrical chi:15; 2) for romanova2018 par1:0.15 par2:0.92" )
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
				BasicDiskBinaryOptions::roptInitializer(vm),
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

