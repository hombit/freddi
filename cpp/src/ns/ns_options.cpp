#include <ns/ns_options.hpp>

NeutronStarOptions::NeutronStarOptions(const po::variables_map &vm):
		NeutronStarArguments(
				vm["nsprop"].as<std::string>(),
				varToOpt<double>(vm, "freqx"),
				varToOpt<double>(vm, "Rx"),
				vm["Bx"].as<double>(),
                vm["Rm_type"].as<std::string>(),             
                vm["h2r_bozzo"].as<double>(),             
                vm["chi_oblique"].as<double>(),
				vm["hotspotarea"].as<double>(),
				vm["epsilonAlfven"].as<double>(),
				vm["inversebeta"].as<double>(),
				vm["Rdead"].as<double>(),
				vm["fptype"].as<std::string>(),
               
				fpparamsInitializer(vm),
				vm["kappattype"].as<std::string>(),
				kappatparamsInitalizer(vm),
				vm["nsgravredshift"].as<std::string>()) {}

				
/*
pard NeutronStarOptions::Rm_typeInitializer(const po::variables_map& vm) {
	const auto Rm_type = vm["Rm_type"].as<std::string>();

	if (Rm_type == "basic") {
		return {};
	}
	if (Rm_type == "bozzo2018") {
		return {};
	}
	throw po::invalid_option_value("Unknown --Rm_type=" + Rm_type);
} 				
*/

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
					"Neutron star properties name: defines geometry (default values of --Rx, --Risco, and --freqx) and accretion->radiation efficiency of NS\n\n"
	 				"Values:\n"
	  				"  dummy: NS accretion->radiation efficiency is R_g * (1 / R_x - 1 / 2R_in), default --freqx is 0, default Rx is 1e6, default Risco is Kerr value\n"
	   				"  newt: NS accretion->radiation efficiency is a function of NS frequency, calculated in Newtonian mechanics (see Lipunova+2021), that's why --freqx must be specified explicitly\n"
					"  sibgatullin-sunyaev2000: NS accretion->radiation efficiency and default values of Rx and Risco are functions of NS frequency, calculated for a specific equation of state for a NS with weak magnetic field (Sibgatullin & Sunyaev, 2000, Astronomy Letters, 26, 699), that's why --freqx must be specified explicitly")
			( "freqx", po::value<double>(), "Accretor rotation frequency, Hz. This parameter is not linked to --kerr, which could be reconciled manually (currently, --kerr is not needed for freddi-ns)" )
			( "Rx", po::value<double>(), "Accretor radius, cm" )
			( "Bx", po::value<double>()->required(), "Accretor polar magnetic induction, G" )
            ( "Rm_type", po::value<std::string>()->default_value(default_Rm_type), 
                    "Definition of magnetosphere radus used in code\n\n"
                    "Values:\n"
                    "  alfven: magnetosphere radius is usual Alfven radius, i.e. Rm = (mu_magn^2 / M_dot / sqrt(G Mx)) ^ 2/7, inclination --chi_oblique not included"
                    "  bozzo: magnetosphere raduis is calculated as described in Bozzo et al. A&A 617, A126 (2018), equation 20, inclination --chi_oblique may be specified" )
            ( "h2r_bozzo", po::value<double>()->default_value(default_h2r_bozzo), "half-thickness to radius relation. Currently used in --bozzo method. " )
            ( "chi_oblique", po::value<double>()->default_value(default_chi_oblique), "Angle between NS magnetic axis and accretion disc axis, deg" )
			( "hotspotarea", po::value<double>()->default_value(default_hotspotarea), "Total area of the region on the accretor radiating beacuse of accretion, normalized by the accretor surface area" )
			( "epsilonAlfven", po::value<double>()->default_value(default_epsilonAlfven), "Magnetosphere radius in units of the Alfven radius, which is defined as (mu^4/G/M/sqrt(Mdot))^(1/7)" ) 
			( "inversebeta", po::value<double>()->default_value(default_inversebeta), "Not currently in use" )
			( "Rdead", po::value<double>()->default_value(default_Rdead), "Maximum inner radius of the disk that can be achieved, cm" )
			( "fptype", po::value<std::string>()->default_value(default_fptype),
					"Scenario to determine the fraction fp of accretied mass. The rest of the disk inner accretion rate is propelled away.\n\n"
					"Values:\n"
	 				"  no-outflow: the matter reaching the inner disk radius always falls onto NS, fp = 1\n"
	  				"  propeller: the matter always flows away, fp = 0\n"
	   				"  corotation-block: like 'no-outflow' when the inner disk radius is smaller than the corotation radius, like 'propeller' otherwise\n"
					"  geometrical: experimental. Generalization of 'corotation-block' for the case of misaligned NS magnetic axis. Requires --fp-geometrical-chi to be specified\n"
					"  eksi-kutlu2010: Under construction\n"
	 				"  romanova2018: fp is an analytical function of the fastness, found from MHD simulations by Romanova et al. (2018, NewA, 62, 94): fp = 1 - par1*fastness^par2. This requires --romanova2018-par1 and --romanova2018-par2 to be specified" )
			( "fp-geometrical-chi", po::value<double>(), "angle between the disk rotation axis and the NS magnetic axis, used for --fptype=geometrical, degrees" )
			( "romanova2018-par1", po::value<double>(), "par1 value for --fptype=romanova2018 and --kappattype=romanova2018" )
			( "romanova2018-par2", po::value<double>(), "par2 value for --fptype=romanova2018 and --kappattype=romanova2018" )
			( "kappattype", po::value<std::string>()->default_value(default_kappat_type),
					"kappa_t describes how strong is the interaction between the NS magnitosphere and disk: total (accelerating) magnetic torque applied to the disc is kappa_t(R) * mu^2 / R^3.\n\n"
					"Values:\n"
					"  const: doesn't depend on radius, kappa_t = value. Requires --kappat-const-value to be specified\n"
					"  corstep: kappa_t can be different inside and outside the corotation radius. Requires --kappat-corstep-in and --kappat-corstep-out to be specified\n"
					"  romanova2018: experimental. Similar to corstep option, but the outside value is reduced by the portion taken away by the outflow (see Table 2 of Romanova+2018, NewA, 62, 94). Requires --kappat-romanova2018-in, --kappat-romanova2018-out --romanova2018-par1 and --romanova-par2 to be specified" )
			( "kappat-const-value", po::value<double>()->default_value(default_kappat_value), "kappa_t value for --kappattype=const" )
			( "kappat-corstep-in", po::value<double>()->default_value(default_kappat_value), "kappa_t value inside the corotation radius for --kappattype=corstep" )
			( "kappat-corstep-out", po::value<double>()->default_value(default_kappat_value), "kappa_t value outside the corotation radius for --kappattype=corstep" )
			( "kappat-romanova2018-in", po::value<double>()->default_value(default_kappat_value), "kappa_t value inside the corotation radius for --kappattype=romanova2018" )
			( "kappat-romanova2018-out", po::value<double>()->default_value(default_kappat_value), "kappa_t value outside the corotation radius for --kappattype=romanova2018" )
			( "nsgravredshift", po::value<std::string>()->default_value(default_ns_grav_redshift),
					"Neutron star gravitational redshift flag.\n\n"
					"Values:\n"
					"  off: gravitational redshift is not taken into account\n"
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


NeutronStarDiskStructureOptions::NeutronStarDiskStructureOptions(const po::variables_map& vm,
																 const NeutronStarArguments& ns_args,
																 const BasicDiskBinaryArguments& bdb_args):
		NeutronStarDiskStructureArguments(
				ns_args,
				bdb_args,
				vm["opacity"].as<std::string>(),
				vm["Mdotout"].as<double>(),
				vm["boundcond"].as<std::string>(),
				vm["Thot"].as<double>(),
				std::pow(vm["Qirr2Qvishot"].as<double>(), 0.25),
				vm["initialcond"].as<std::string>(),
				varToOpt<double>(vm, "F0"),
				varToOpt<double>(vm, "Mdisk0"),
				varToOpt<double>(vm, "Mdot0"),
				varToOpt<double>(vm, "powerorder"),
				varToOpt<double>(vm, "gaussmu"),
				varToOpt<double>(vm, "gausssigma"),
				vm["windtype"].as<std::string>(),
				DiskStructureOptions::windparamsInitializer(vm)) {}

po::options_description NeutronStarDiskStructureOptions::description() {
	auto non_ns_od = DiskStructureOptions::description();
	po::options_description od(DiskStructureOptions::caption);
	for (const auto &non_ns_option : non_ns_od.options()) {
		if (non_ns_option->long_name() == "initialcond") {
			od.add_options()
					( "initialcond", po::value<std::string>()->default_value(default_initialcond), (non_ns_option->description() + "  quasistat-ns: Distibution of the initial viscous torque in the disc is  F0 = Mdot0 * (h_out - h_in) / h_out * h_in / oprel.f_F(h_in / h_out), where f_F is taken from Lipunova & Shakura (2000)\n").c_str() )
					;
		} else {
			od.add(non_ns_option);
		}
	}
	return od;
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
			( "angulardistns", po::value<std::string>()->default_value(default_angular_dist_ns), "Flag to calculate angular distribution the NS emission. Values: isotropic, plane" )
			;
	return od;
}


FreddiNeutronStarOptions::FreddiNeutronStarOptions(const po::variables_map &vm) {
	ns.reset(new NeutronStarOptions(vm));

	general.reset(new GeneralOptions(vm));
	basic.reset(new NeutronStarBasicDiskBinaryOptions(vm, *ns));
	disk.reset(new NeutronStarDiskStructureOptions(vm, *ns, *basic));
	irr_ns.reset(new NeutronStarSelfIrradiationOptions(vm, *disk));
	irr = irr_ns;
	flux.reset(new FluxOptions(vm));
	calc.reset(new CalculationOptions(vm));
}

po::options_description FreddiNeutronStarOptions::description() {
	po::options_description desc("Freddi NS: numerical calculation of accretion disk evolution");
	desc.add(GeneralOptions::description());
	desc.add(NeutronStarBasicDiskBinaryOptions::description());
	desc.add(NeutronStarDiskStructureOptions::description());
	desc.add(NeutronStarOptions::description());
	desc.add(NeutronStarSelfIrradiationOptions::description());
	desc.add(FluxOptions::description());
	desc.add(CalculationOptions::description());
	return desc;
}
