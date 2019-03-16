#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/arrayobject.h>

#include <util.hpp>

#include "converters.hpp"


using namespace boost::python;
namespace np = boost::python::numpy;


struct VectorToListConverter {
	static PyObject* convert(const vecd& v) {
		list l;
		for (auto &x : v) {
			l.append(x);
		}
		return incref(l.ptr());
	}
};

template<typename T>
struct VectorToNumpyConverter {
	static PyObject* convert(const std::vector<T> v) {
		np::ndarray a = np::from_data(
				v.data(), np::dtype::get_builtin<T>(),
				make_tuple(v.size()), make_tuple(sizeof(T)),
				object());
		return incref(a.copy().ptr());
	}
};

template<typename T>
struct NumpyToVectorConverter {
	NumpyToVectorConverter() {
		converter::registry::push_back(convertible, construct, type_id<std::vector<T> >());
	}
	static void* convertible(PyObject* obj) {
		return PyArray_Check(obj) ? obj : nullptr;
	}
	static void construct(PyObject* obj, converter::rvalue_from_python_stage1_data* data) {
		auto arr = reinterpret_cast<PyArrayObject*>(obj);
		const auto size = PyArray_SIZE(arr);
		T* arr_data = reinterpret_cast<T*>(PyArray_DATA(arr));

		auto storage = reinterpret_cast<converter::rvalue_from_python_storage<std::vector<T> >*>(data)->storage.bytes;
		new(storage) std::vector<T>(arr_data, arr_data + size);
		data->convertible = storage;
	}
};


void register_converters() {
	to_python_converter<vecd, VectorToNumpyConverter<double>, false>();
//	NumpyToVectorConverter<double>();
}
