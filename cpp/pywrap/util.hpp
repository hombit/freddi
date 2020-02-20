#ifndef FREDDI_PYWRAP_UTIL_HPP
#define FREDDI_PYWRAP_UTIL_HPP



#include <optional>

#include <boost/python.hpp>

using namespace boost::python;


template<class T>
std::optional<T> objToOpt(const object& obj) {
	static const auto None = object().ptr();
	if (obj.ptr() == None) {
		return {};
	}
	return extract<T>(obj);
}


#endif //FREDDI_PYWRAP_UTIL_HPP
