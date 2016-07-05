#include "orbit.hpp"


double r_out_func(const double Mx, const double Mopt, const double P){
	const double semiAxis = cbrt( GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * ( Mx + Mopt ) * P * P / ( 4. * M_PI*M_PI ) );	
	const double q = cbrt(Mx / Mopt);
	const double roche = semiAxis * 0.49 * q*q / ( 0.6 * q*q + log(1.+q) ); // Volume radius of Roche lobe, approximation from Eggleton, P. P. 1983, ApJ, 268, 368
	return 0.8 * roche;
}

double r_in_func(double Mx, double A){ // From «Black Hole Accretion Disks», A.44 (p. 530)
			double Z1 = 1. + cbrt( (1.-A*A) ) * ( cbrt( (1.+A) ) + cbrt( (1.-A) ) );
			double Z2 = sqrt( 3.*A*A + Z1*Z1 );
			return GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT) * ( 3. + Z2 - sqrt( (3.-Z1) * (3.+Z1+2.*Z2) ) );
};
