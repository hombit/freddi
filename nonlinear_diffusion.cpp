#include "nonlinear_diffusion.hpp"



double mean_square_rel(const vecd &A, const vecd &B, int first, int last){
	double rv = 0;
	double odds;
	for ( int i = first; i <= last; ++i ){
		odds = ( A.at(i) - B.at(i) ) / A.at(i);
		rv += odds*odds;
	}
	return sqrt(rv) / (last-first+1);
}


double max_dif_rel(const vecd &A, const vecd &B, int first, int last){
	double max = 0.;
	for ( int i = first; i <= last; ++i ){
		double x = fabs ( ( A.at(i) - B.at(i) ) / A.at(i) );
		if ( x > max )
			max = x;
	}
	return max;
}



// \frac{dw}{dt}=\frac{d^2y}{dx^2}, y=y(x,t) — ?, w = w (x,y)

void nonlenear_diffusion_nonuniform_1_2 (const double tau,
										 const double eps, // reletive error for w
										 const double left_bounder_cond, // y(left_border,Time+tau) = left_bounder_cond
										 const double right_bounder_cond, // \frac{y(right_border,Time+tau)}{dx} = right_bounder_cond
										 std::function<vecd (const vecd &, const vecd &, unsigned int, unsigned int)> wunc, // first argument is array of x_i, second — array of y(x_i,t); return value — array of w(x_i,y_i)
										 const vecd &x, // array with (non)uniform grid
										 vecd &y// array with initial coundition and for results
									 	){
	vecd y_init(y);
	const int N = fmin(x.size(), y.size()); // N_x+1
	const auto W = wunc(x, y, 1, N-1);
	vecd K_0(N), K_1(N), CC(N), frac(N), a(N), b(N), f(N);
	for ( int i = 1; i < N-1; ++i ){
		a.at(i) = 2. * ( x.at(i+1) - x.at(i) ) / ( x.at(i+1) - x.at(i-1) );
		b.at(i) = 2. * ( x.at(i) - x.at(i-1) ) / ( x.at(i+1) - x.at(i-1) );
		frac.at(i) = ( x.at(i+1) - x.at(i) ) * ( x.at(i) - x.at(i-1) ) / tau;
	}
    frac.at(N-1) = ( x.at(N-1) - x.at(N-2) ) * ( x.at(N-1) - x.at(N-2) ) * 0.5 / tau;
    for ( int i = 1; i < N; ++i ){
        f.at(i) = frac.at(i) * W.at(i);
		K_1.at(i) = frac.at(i) * W.at(i) / y.at(i);
		CC.at(i) = K_0.at(i) = K_1.at(i) * (1. + 2.*eps);
    }
	auto iteration = [&](vecd &K) -> void{ // [&] <-> [&wunc, &x, &y, &frac, &a, &b, &f, N, left_bounder_cond, right_bounder_cond]
		vecd alpha(N), beta(N);
		alpha.at(1) = 0.;
		beta.at(1) = left_bounder_cond;
		for ( int i = 1; i < N-1; ++i ){
			double c = 2. + K.at(i);
			alpha.at(i+1) = b.at(i) / ( c - alpha.at(i) * a.at(i) );
			beta.at(i+1) = ( beta.at(i) * a.at(i) + f.at(i) ) / ( c - alpha.at(i) * a.at(i) );
		}
		// y.at(N-1) = ( (x.at(N-1) - x.at(N-2) ) * right_bounder_cond + beta.at(N-1) ) / ( 1. - alpha.at(N-1) );
		y.at(N-1) = ( (x.at(N-1) - x.at(N-2) ) * right_bounder_cond + f.at(N-1) + beta.at(N-1) ) / ( 1. + K.at(N-1) - alpha.at(N-1) );
		for ( int i = N-2; i > 0; --i )
			y.at(i) = alpha.at(i+1) * y.at(i+1) + beta.at(i+1);
		y.at(0) = left_bounder_cond;
		auto WW = wunc(x, y, 1, N-1);
		for ( int i = 1; i < N-1; ++i )
			K.at(i) = frac.at(i) * WW.at(i) / y.at(i);
	};

	bool flag = false;	int j = 0;	double delta;	vecd D_1(N), D_2(N);
	while( max_dif_rel(K_1, K_0, 1, N-2) > eps ){
		if ( max_dif_rel(K_1, CC, 1, N-2) > 0. and flag == false ){
			K_0 = K_1;
			iteration(K_1);

			j++;
			if ( j % 2 == 0 )
				CC = K_0;
			if ( j % 4 == 1 )
				delta = max_dif_rel (K_1, K_0, 1, N-2);
			if ( j % 4 == 3 and max_dif_rel (K_1, K_0, 1, N-2) >= delta ){
				flag = true;
			}
		} else{
			throw std::runtime_error("Divergence in nonlinear_diffusion");
		}
	}
}



void nonlenear_diffusion_nonuniform_2_2 (const double tau,
										 const double eps, // reletive error for w
										 const double left_bounder_cond, // // \frac{y(left_border,Time+tau)}{dx} = left_bounder_cond
										 const double right_bounder_cond, // \frac{y(right_border,Time+tau)}{dx} = right_bounder_cond
										 std::function<vecd (const vecd &, const vecd &, unsigned int, unsigned int)> wunc, // first argument is array of x_i, second — array of y(x_i,t); return value — array of w(x_i,y_i)
										 const vecd &x, // array with (non)uniform grid
										 vecd &y// array with initial coundition and for results
									 	){
	vecd y_init(y);
	const int N = fmin(x.size(), y.size()); // N_x+1
	const auto W = wunc(x, y, 1, N-1);
	vecd K_0(N), K_1(N), CC(N), frac(N), a(N), b(N), f(N);
	for ( int i = 1; i < N-1; ++i ){
		frac.at(i) = ( x.at(i+1) - x.at(i) ) * ( x.at(i) - x.at(i-1) ) / tau;
		a.at(i) = 2. * ( x.at(i+1) - x.at(i) ) / ( x.at(i+1) - x.at(i-1) );
		b.at(i) = 2. * ( x.at(i) - x.at(i-1) ) / ( x.at(i+1) - x.at(i-1) );
		f.at(i) = frac.at(i) * W.at(i);
		K_1.at(i) = frac.at(i) * W.at(i) / y.at(i);
		CC.at(i) = K_0.at(i) = K_1.at(i) * 2. + 10.*eps;
	}
	auto iteration = [&](vecd &K) -> void{ // [&] <-> [&wunc, &x, &y, &frac, &a, &b, &f, N, left_bounder_cond, right_bounder_cond]
		vecd alpha(N), beta(N);
		alpha.at(1) = 1.;
		beta.at(1) = - (x.at(1) - x.at(0)) * left_bounder_cond;
		for ( int i = 1; i < N-1; ++i ){
			double c = 2. + K.at(i);
			alpha.at(i+1) = b.at(i) / ( c - alpha.at(i) * a.at(i) );
			beta.at(i+1) = ( beta.at(i) * a.at(i) + f.at(i) ) / ( c - alpha.at(i) * a.at(i) );
		}
		y.at(N-1) = ( (x.at(N-1) - x.at(N-2) ) * right_bounder_cond + beta.at(N-1) ) / ( 1. - alpha.at(N-1) );
		for ( int i = N-2; i > 0; --i )
			y.at(i) = alpha.at(i+1) * y.at(i+1) + beta.at(i+1);
		y.at(0) = left_bounder_cond;
		auto WW = wunc(x, y, 1, N-1);
		for ( int i = 1; i < N-1; ++i )
			K.at(i) = frac.at(i) * WW.at(i) / y.at(i);
	};

	bool flag = false;	int j = 0;	double delta;	vecd D_1(N), D_2(N);
	while( max_dif_rel(K_1, K_0, 1, N-2) > eps ){
		if ( max_dif_rel(K_1, CC, 1, N-2) > 0. and flag == false ){
			K_0 = K_1;
			iteration(K_1);

			j++;
			if ( j % 2 == 0 )
				CC = K_0;
			if ( j % 4 == 1 )
				delta = max_dif_rel (K_1, K_0, 1, N-2);
			if ( j % 4 == 3 and max_dif_rel (K_1, K_0, 1, N-2) >= delta ){
				flag = true;
			}
		} else{
			throw std::runtime_error("Divergence in nonlinear_diffusion");
		}
	}
}
