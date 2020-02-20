#ifndef FREDDI_PYTHON_CONVERTERS_HPP
#define FREDDI_PYTHON_CONVERTERS_HPP

#include <boost/python.hpp>

#include <util.hpp>


using object = boost::python::object;

pard mapping_to_map(const object& obj);

void register_converters();


#endif //FREDDI_PYTHON_CONVERTERS_HPP
