#include "opacity_related.hpp"


OpacityRelated::OpacityRelated(
	const std::string &opacity_type,
	double Mx,
	double alpha,
	double mu
):
	type(opacity_type),
	Mx(Mx),
	alpha(alpha),
	mu(mu)
{
	GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx;

	if ( type == "Kramers" ){
		init_Kramers();
	} else if ( type == "OPAL" ){
		init_OPAL();
	} else{
		throw std::invalid_argument(opacity_type);
	}
}


void OpacityRelated::init_Kramers(){
	m = 0.3;
	n = 0.8;
	varkappa0 = 5e24;

	// tau_0 = 1e3
	Pi1 = 6.31217;
	Pi2 = 0.51523;
	Pi3 = 1.13050;
	Pi4 = 0.39842;

	Pi_Sigma = pow(Pi1, 0.05) * pow(Pi2, 0.1) * pow(Pi3, 0.8) * pow(Pi4, 0.1);
	Pi_Height = pow(Pi1, 19./40.) * pow(Pi2, -0.05) * pow(Pi3, 0.1) * pow(Pi4, -0.05);

	D = 2.41e34 * pow(alpha, 0.8) * (Mx/GSL_CONST_CGSM_SOLAR_MASS) * pow(mu/0.5, -0.75) / Pi_Sigma * pow(varkappa0/1e22, 0.1);

	Height_exp_F = 3./20.;
	Height_exp_R = 1./8.;
	Height_coef = 0.020 * pow(1e17, -Height_exp_F) * pow(GM, -Height_exp_F/2.) * pow(1e10, -Height_exp_F/2.) * pow(alpha, -0.1) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -3./8.) * pow(mu/0.6, -3./8.) * Pi_Height * pow(varkappa0/5e24, 0.05);


	a0 = 1.430;
	a1 = -0.46;
	a2 = 0.03;
	k = 3.5;
	l = 6.0;
}


void OpacityRelated::init_OPAL(){
	m = 1./3.;
	n = 1.;
	varkappa0 = 1.5e20;

	// tau_0 = 1e3
	Pi1 = 5.53450;
	Pi2 = 0.55220;
	Pi3 = 1.14712;
	Pi4 = 0.44777;

	Pi_Sigma = pow(Pi1, 1./18.) * pow(Pi2, 1./9.) * pow(Pi3, 7./9.) * pow(Pi4, 1./9.);
	Pi_Height = pow(Pi1, 17./36.) * pow(Pi2, -1./18.) * pow(Pi3, 1./9.) * pow(Pi4, -1./18.);

	D = 2.12e37 * pow(alpha, 7./9.) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, 10./9.) * pow(mu/0.5, -13./18.) / Pi_Sigma * pow(varkappa0/1e22, 1./9.);
	
	Height_exp_F = 1./6.;
	Height_exp_R = 1./12.;
	Height_coef = 0.021 * pow(1e17, -Height_exp_F) * pow(GM, -Height_exp_F/2.) * pow(1e10, -Height_exp_F/2.) * pow(alpha, -1./9.) * pow(Mx/GSL_CONST_CGSM_SOLAR_MASS, -13./36.) * pow(mu/0.6, -13./36.) * Pi_Height * pow(varkappa0/1.5e20, 1./18.);

	
	a0 = 1.3998;
	a1 = -0.425;
	a2 = 0.0249;
	k = 11./3.;
	l = 19./3.;
}


double OpacityRelated::Height(double R, double F) const {
	return R * Height_coef * pow(F, Height_exp_F) * pow(R/1e10, Height_exp_R - Height_exp_F/2.);
}


double OpacityRelated::f_F(double xi) const {
	return a0 * xi + a1 * pow(xi, k) + a2 * pow(xi, l);
}
