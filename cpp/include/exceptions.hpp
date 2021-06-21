#ifndef FREDDI_EXCEPTIONS_HPP
#define FREDDI_EXCEPTIONS_HPP

#include <exception>

class RadiusCollapseException: public std::exception {
public:
	virtual const char* what() const noexcept override {
		return "Rout <= Rin";
	}
};

#endif //FREDDI_EXCEPTIONS_HPP
