#ifndef _OPACITY_RELATED_HPP
#define _OPACITY_RELATED_HPP


#include <algorithm> // std::none_of
#include <array>
#include <cmath>
#include <stdexcept> // std::invalid_argument
#include <string>

#include "gsl_const_cgsm.h"


class OpacityRelated{
private:
	void init_Kramers();
	void init_OPAL();

public:
	OpacityRelated(
		const std::string &opacity_type,
		double Mx,
		double alpha,
		double mu
	);
	~OpacityRelated(){};

	const std::array<std::string, 2> supported_types {{ "Kramers", "OPAL" }};
	const std::string type;

	const double Mx, alpha, mu;
	double GM;
	double m, n, varkappa0, Pi1, Pi2, Pi3, Pi4, Pi_Sigma, Pi_Height, D, Height_exp_F, Height_exp_R, Height_coef;
	double a0, a1, a2, k, l;

	double Height(double R, double F) const;
	double f_F(double xi) const;
};


#endif // _OPACITY_RELATED_HPP
