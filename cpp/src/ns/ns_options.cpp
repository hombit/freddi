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
				fpparamsInitializer(vm),
				vm["kappattype"].as<std::string>(),
				kappatparamsInitalizer(vm),
				vm["nsgravredshift"].as<std::string>()) {}

pard NeutronStarOptions::fpparamsInitializer(const po::variables_map& vm) {
	const auto fptype = vm["fptype"].as<std::string>();

	if (fptype == "no-outflow") {
		return {};
	}
	if (fptype == "propeller") {
		return {};
	}
	if (fptype == "corotation-block") {
		return {};
	}
	if (fptype == "geometrical") {
		if (vm.count("fp-geometrical-chi") == 0) {
			throw po::error("--fp-geometrical-chi is required if --fptype=geometrical");
		}
		return {{"chi", vm["fp-geometrical-chi"].as<double>()}};
	}
	if (fptype == "eksi-kultu2010") {
		return {};
	}
	if (fptype == "romanova2018") {
		if (vm.count("romanova2018-par1") == 0) {
			throw po::error("--romanova2018-par1 is required if --fptype=romanova2018");
		}
		const double par1 = vm["romanova2018-par1"].as<double>();
		if (vm.count("romanova2018-par2") == 0) {
			throw po::error("--romanova2018-par2 is required if --fptype=romanova2018");
		}
		const double par2 = vm["romanova2018-par2"].as<double>();
		return {
			{"par1", par1},
			{"par2", par2},
		};
	}

	throw po::invalid_option_value("Unknown --fptype=" + fptype);
}

pard NeutronStarOptions::kappatparamsInitalizer(const po::variables_map &vm) {
	const auto kappattype = vm["kappattype"].as<std::string>();

	if (kappattype == "const") {
		if (vm.count("kappat-const-value") == 0) {
			throw po::error("--kappat-const-value is required if --kappattype=const");
		}
		return {{"value", vm["kappat-const-value"].as<double>()}};
	}
	if (kappattype == "corstep") {
		if (vm.count("kappat-corstep-in") == 0) {
			throw po::error("--kappat-corstep-in is required if --kappattype=corstep");
		}
		const double in = vm["kappat-corstep-in"].as<double>();
		if (vm.count("kappat-corstep-out") == 0) {
			throw po::error("--kappat-corstep-out is required if --kappattype=corstep");
		}
		const double out = vm["kappat-corstep-out"].as<double>();
		return {
			{"in", in},
			{"out", out},
		};
	}
	if (kappattype == "romanova2018") {
		if (vm.count("kappat-romanova2018-in") == 0) {
			throw po::error("--kappat-romanova2018-in is required if --kappattype=romanova2018");
		}
		const double in = vm["kappat-romanova2018-in"].as<double>();
		if (vm.count("kappat-romanova2018-out") == 0) {
			throw po::error("--kappat-romanova2018-out is required if --kappattype=romanova2018");
		}
		const double out = vm["kappat-romanova2018-out"].as<double>();
		if (vm.count("romanova2018-par1") == 0) {
			throw po::error("--romanova2018-par1 is required if --kappattype=romanova2018");
		}
		const double par1 = vm["romanova2018-par1"].as<double>();
		if (vm.count("romanova2018-par2") == 0) {
			throw po::error("--romanova2018-par2 is required if --kappattype=romanova2018");
		}
		const double par2 = vm["romanova2018-par2"].as<double>();
		return {
			{"in", in},
			{"out", out},
			{"par1", par1},
			{"par2", par2},
		};
	}

	throw po::invalid_option_value("Unknown --kappattype=" + kappattype);
}

po::options_description NeutronStarOptions::description() {
	po::options_description od("Parameters of accreting neutron star");
	od.add_options()
			( "nsprop", po::value<std::string>()->default_value(default_nsprop),
					"Neutron star geometry and radiation properties name and specifies default values of --Rx, --Risco and --freqx\n\n"
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
			( "fptype", po::value<std::string>()->default_value(default_fptype),
					"??? Accretor Mdot fraction mode fp.\n\n"
					"Values:\n"
	 				"  no-outflow: all the matter passing inner disk radius falling onto neutron star, fp = 1\n"
	  				"  propeller: all the matter flows away from both disk and neutron star, fp = 0\n"
	   				"  corotation-block: like 'no-otflow' when Alfven radius is smaller than corotation radius, like 'propeller' otherwise\n"
					"  geometrical: generalisation of 'corotation-block' for the case of not co-directional of disk rotation axis and neutron star magnetic field axis. Requires --fp-geometrical-chi to be specified\n"
					"  eksi-kutlu2010: ???\n"
	 				"  romanova2018: ???, requires --romanova2018-par1 and --romanova2018-par2 to be specified" )
			( "fp-geometrical-chi", po::value<double>(), "angle between disk rotation axis and neutron star magnetic axis for --fptype=geometrical, degrees" )
			( "romanova2018-par1", po::value<double>(), "??? par1 value for --fptype=romanova2018 and --kappattype=romanova2018" )
			( "romanova2018-par2", po::value<double>(), "??? par2 value for --fptype=romanova2018 and --kappattype=romanova2018" )
			( "kappattype", po::value<std::string>()->default_value(default_kappat_type),
					"kappa_t describes how strong is interaction between neutron star magnitosphere and disk, magnetic torque is kappa_t(R) * mu^2 / R^3. This parameter describes type of radial destribution of this parameter\n\n"
					"Values:\n"
					"  const: doesn't depend on radius, kappa_t = value. Requires --kappat-const-value to be specified\n"
					"  corstep: kappa_t is 'in' inside corotation radius, and 'out' outside. Requires --kappat-corstep-in and --kappat-corstep-out to be specified\n"
					"  romanova2018: similar to corstep option, but the outside value is reduced by the portion taken away by wind (see Table 2 of Romanova+2018,NewA,62,94). Requires --kappat-romanova2018-in, --kappat-romanova2018-out --romanova2018-par1 and --romanova-par2 to be specified" )
			( "kappat-const-value", po::value<double>()->default_value(default_kappat_value), "kappa_t value for --kappattype=const" )
			( "kappat-corstep-in", po::value<double>()->default_value(default_kappat_value), "kappa_t value inside corotation radius for --kappattype=corstep" )
			( "kappat-corstep-out", po::value<double>()->default_value(default_kappat_value), "kappa_t value outside corotation radius for --kappattype=corstep" )
			( "kappat-romanova2018-in", po::value<double>()->default_value(default_kappat_value), "kappa_t value inside corotation radius for --kappattype=romanova2018" )
			( "kappat-romanova2018-out", po::value<double>()->default_value(default_kappat_value), "kappa_t value outside corotation radius for --kappattype=romanova2018" )
			( "nsgravredshift", po::value<std::string>()->default_value(default_ns_grav_redshift),
					"Neutron star gravitational redshift type.\n\n"
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
				varToOpt<double>(vm, "alphacold"),
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
				vm["h2rcold"].as<double>(),
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
