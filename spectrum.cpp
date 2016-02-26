#include "spectrum.hpp"

double Luminosity( const std::vector<double> &R, const std::vector<double> &T, double min_nu, double max_nu, int Nnu ){
	const int NR = fmin(R.size(), T.size());
	const double step_nu = Nnu > 1.  ?  ( max_nu - min_nu ) / (Nnu-1.)  :  1.;
	double L = 0;
	for ( int i_nu = 0; i_nu < Nnu; ++i_nu ){
		const double nu = min_nu + step_nu * i_nu;
		double Inu = 0;
		for ( int i_R = 0; i_R < NR; ++i_R ){
			double stepR;
			if ( i_R == 0               ) stepR = R.at(i_R+1) - R.at(i_R  );
			if ( i_R == NR-1            ) stepR = R.at(i_R  ) - R.at(i_R-1);
			if ( i_R > 1 and i_R < NR-1 ) stepR = R.at(i_R+1) - R.at(i_R-1);
			const double Bnu = 2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H * nu * nu * nu / GSL_CONST_CGSM_SPEED_OF_LIGHT / GSL_CONST_CGSM_SPEED_OF_LIGHT / ( exp( nu*GSL_CONST_CGSM_PLANCKS_CONSTANT_H / GSL_CONST_CGSM_BOLTZMANN / T.at(i_R) ) - 1. );
			Inu += .5 * Bnu * 2. * M_PI * R.at(i_R) * stepR;
		}
		if ( (i_nu == 0 or i_nu == Nnu-1) and Nnu > 1. ){
			L += Inu / 2.;
		} else{
			L += Inu;
		}
	}
	L *= 2. * M_PI * step_nu;
	return L;
}



double I_lambda( const std::vector<double> &R, const std::vector<double> &T, double lambda ){
	double I = 0;
	int NR = fmin(R.size(), T.size());
	for ( int i_R = 1; i_R < NR; ++i_R ){
		double stepR;
		if ( i_R == 0              ) stepR = R.at(i_R+1) - R.at(i_R  );
		if ( i_R == NR-1           ) stepR = R.at(i_R  ) - R.at(i_R-1);
		if ( i_R > 1 and i_R < NR-1 ) stepR = R.at(i_R+1) - R.at(i_R-1);
		const double B_lambda =  2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT / pow(lambda,5.) / ( exp( GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_PLANCKS_CONSTANT_H / lambda / GSL_CONST_CGSM_BOLTZMANN / T.at(i_R) ) - 1. );
		I += .5 * B_lambda * 2. * M_PI * R.at(i_R) * stepR;
	}
	return I;
}


// Code by Galina Lipunova:
/* General Relativity effects are included in the structure of the disk
   (Page & Thorne 1974; Riffert & Herold 1995). metric = "GR"
*/
double T_GR( const double r1, const double ak, const double Mx, const double Mdot, const double r_in  ){
    const double GM = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * Mx;
    const double rg = GM  / (GSL_CONST_CGSM_SPEED_OF_LIGHT*GSL_CONST_CGSM_SPEED_OF_LIGHT);
    const double x = sqrt(r1 / rg);
    const double x0 = sqrt(r_in/rg);
    
    const double x1 = 2. * cos ((acos(ak)-M_PI)/3.);
    const double x2 = 2. * cos ((acos(ak)+M_PI)/3.);
    const double x3 = -2. * cos (acos(ak)/3.);
    const double a = 3. * (x1-ak)*(x1-ak) * log((x-x1)/(x0-x1))/x1/(x1-x2)/(x1-x3);
    const double b = 3. * (x2-ak)*(x2-ak) * log((x-x2)/(x0-x2))/x2/(x2-x1)/(x2-x3);
    const double c = 3. * (x3-ak)*(x3-ak) * log((x-x3)/(x0-x3))/x3/(x3-x1)/(x3-x2);
   
    return( pow(
        (3.*Mdot * pow(GSL_CONST_CGSM_SPEED_OF_LIGHT,6.) / (8.*M_PI*GM*GM)) *
            (x-x0-1.5*ak*log(x/x0)-a-b-c) / ( pow(x,4.)*(x*x*x-3.*x+2.*ak) ),
        0.25)
    );
}
