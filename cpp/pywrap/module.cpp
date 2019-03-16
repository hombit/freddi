#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "arguments.hpp"
#include "converters.hpp"
#include "freddi_evolution.hpp"
#include "freddi_state.hpp"


namespace np = boost::python::numpy;


BOOST_PYTHON_MODULE(_freddi) {
	Py_Initialize();
	np::initialize();

	register_converters();

	wrap_arguments();
	wrap_state();
	wrap_evolution();
}
