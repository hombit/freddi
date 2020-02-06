#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "converters.hpp"
#include "pywrap_arguments.hpp"
#include "pywrap_freddi_evolution.hpp"
#include "pywrap_freddi_state.hpp"


namespace np = boost::python::numpy;


BOOST_PYTHON_MODULE(_freddi) {
	Py_Initialize();
	np::initialize();

	register_converters();

	wrap_arguments();
	wrap_state();
	wrap_evolution();
}
