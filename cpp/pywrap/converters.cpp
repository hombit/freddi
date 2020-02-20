#include <map>
#include <string>
#include <vector>

#include <boost/python/numpy.hpp>
#include <numpy/arrayobject.h>

#include "converters.hpp"


using namespace boost::python;
namespace np = boost::python::numpy;


pard mapping_to_map(const object& obj) {
	if (!PyMapping_Check(obj.ptr())) {
		PyErr_SetString(PyExc_TypeError, "windparams argument must be a mapping");
		throw error_already_set();
	}
	pard m;
	stl_input_iterator<object> begin(obj), end;
	for (auto key = begin; key != end; ++key) {
		m[extract<std::string>(*key)] = extract<double>(obj[*key]);
	}
	return m;
}


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

template<typename CTYPE, int NPTYPE>
struct NumpyToVectorConverter {
	NumpyToVectorConverter() {
		converter::registry::push_back(convertible, construct, type_id<std::vector<CTYPE> >());
	}
	static void* convertible(PyObject* obj) {
		if (!PyArray_Check(obj)) {
			return nullptr;
		}
		const auto arr = reinterpret_cast<PyArrayObject*>(obj);
		const auto dtype = PyArray_DTYPE(arr);
		if (dtype->type_num != NPTYPE) {
			return nullptr;
		}
		return obj;
	}
	static void construct(PyObject* obj, converter::rvalue_from_python_stage1_data* data) {
		const auto arr = reinterpret_cast<PyArrayObject*>(obj);
		const auto size = PyArray_SIZE(arr);
		CTYPE* arr_data = reinterpret_cast<CTYPE*>(PyArray_DATA(arr));

		auto storage = reinterpret_cast<converter::rvalue_from_python_storage<std::vector<CTYPE> >*>(data)->storage.bytes;
		new(storage) std::vector<CTYPE>(arr_data, arr_data + size);
		data->convertible = storage;
	}
};


template<typename Key, typename T>
struct MapToDictConverter {
	static PyObject* convert(const std::map<Key, T> m) {
		dict d;
		for (const auto &x: m) {
			d[x.first] = x.second;
		}
		return incref(d.ptr());
	}
};


template<typename Key, typename T>
struct DictToMapConverter {
	DictToMapConverter() {
		converter::registry::push_back(convertible, construct, type_id<std::map<Key, T> >());
	}
	static void* convertible(PyObject* obj) {
		if (!PyMapping_Check(obj)) {
			return nullptr;
		}
		return obj;
	}
	static void construct(PyObject* obj, converter::rvalue_from_python_stage1_data* data) {
		handle<> handle(obj);
		dict d(handle);
		std::map<Key, T> m;

		const stl_input_iterator<object> begin(d), end;
		for (auto key = begin; key != end; ++key) {
			extract<Key> get_key(*key);
			if (!get_key.check()) {
				PyErr_SetString(PyExc_TypeError, "Wrong mapping key type");
				throw error_already_set();
			}
			extract<T> get_value(d[*key]);
			if (!get_value.check()) {
				PyErr_SetString(PyExc_TypeError, "Wrong mapping value type");
				throw error_already_set();
			}
			m[get_key()] = get_value();
		}

		auto storage = reinterpret_cast<converter::rvalue_from_python_storage<std::map<Key, T> >*>(data)->storage.bytes;
		new(storage) std::map<Key, T>(std::move(m));
		data->convertible = storage;
	}
};


void register_converters() {
	to_python_converter<vecd, VectorToNumpyConverter<double>, false>();
	NumpyToVectorConverter<double, NPY_DOUBLE>();

	to_python_converter<std::map<std::string, double>, MapToDictConverter<std::string, double>, false>();
	DictToMapConverter<std::string, double>();
}
