#include <vector>

#include <boost/python.hpp>

#include <freddi_state.hpp>
#include <util.hpp>

#include "pywrap_freddi_state.hpp"

using namespace boost::python;


void wrap_state() {
	double (FreddiState::*flux_hot)(double) = &FreddiState::flux;
	double (FreddiState::*flux_cold)(double) = &FreddiState::flux_region<FreddiState::ColdRegion>;
	double (FreddiState::*flux_star)(double, double) = &FreddiState::flux_star;

	class_<FreddiState>("_State", no_init)
	    .add_property("GM", &FreddiState::GM)
	    .add_property("eta", &FreddiState::eta)
	    .add_property("distance", &FreddiState::distance)
		.add_property("Mdot_in", &FreddiState::Mdot_in)
		.add_property("Mdot", &FreddiState::Mdot_in)
		.add_property("Mdot_out", &FreddiState::Mdot_out)
		//.add_property("Lx", &FreddiState::Lx)
        .add_property("Fx", &FreddiState::Fx)
		.add_property("t", &FreddiState::t)
		.add_property("i_t", &FreddiState::i_t)
		.add_property("Nt", &FreddiState::Nt)
		.add_property("Nx", &FreddiState::Nx)
		.add_property("first", &FreddiState::first)
		.add_property("last", &FreddiState::last)
		.add_property("Mdisk", &FreddiState::Mdisk)
		.add_property("Mdot_wind", &FreddiState::Mdot_wind)
		.add_property("h", make_function(&FreddiState::h, return_value_policy<copy_const_reference>()))
		.add_property("R", make_function(&FreddiState::R, return_value_policy<copy_const_reference>()))
		.add_property("F", make_function(&FreddiState::F, return_value_policy<copy_const_reference>()))
		.add_property("W", make_function(&FreddiState::W, return_value_policy<copy_const_reference>()))
		.add_property("Tph", make_function(&FreddiState::Tph, return_value_policy<copy_const_reference>()))
		.add_property("Tph_vis", make_function(&FreddiState::Tph_vis, return_value_policy<copy_const_reference>()))
		.add_property("Tirr", make_function(&FreddiState::Tirr, return_value_policy<copy_const_reference>()))
		.add_property("Kirr", make_function(&FreddiState::Kirr, return_value_policy<copy_const_reference>()))
		.add_property("Sigma", make_function(&FreddiState::Sigma, return_value_policy<copy_const_reference>()))
		.add_property("Height", make_function(&FreddiState::Height, return_value_policy<copy_const_reference>()))
		.add_property("windA", make_function(&FreddiState::windA))
		.add_property("windB", make_function(&FreddiState::windB))
		.add_property("windC", make_function(&FreddiState::windC))
		.add_property("lambdas", make_function(&FreddiState::lambdas, return_value_policy<copy_const_reference>()))
		.def("_flux_hot", flux_hot)
		.def("_flux_cold", flux_cold)
		.def("_flux_star", flux_star)
	;
}
