#include "nonlinear_diffusion.hpp"



double mean_square_rel(const vecd &A, const vecd &B, size_t first, size_t last){
	double rv = 0;
	double odds;
	for ( size_t i = first; i <= last; ++i ){
		odds = ( A[i] - B[i] ) / A[i];
		rv += odds*odds;
	}
	return std::sqrt(rv) / (last-first+1);
}


double max_dif_rel(const vecd &A, const vecd &B, size_t first, size_t last){
	double max = 0.;
	double x;
	for ( size_t i = first; i <= last; ++i ){
		x = std::abs ( ( A[i] - B[i] ) / A[i] );
		if ( x > max ) {
			max = x;
		}
	}
	return max;
}



// \frac{dw}{dt}=\frac{d^2y}{dx^2}, y=y(x,t) — ?, w = w (x,y)

void nonlenear_diffusion_nonuniform_1_2 (const double tau,
										 const double eps, // relative error for w
										 const double left_bounder_cond, // y(left_border,Time+tau) = left_bounder_cond
										 const double right_bounder_cond, // \frac{y(right_border,Time+tau)}{dx} = right_bounder_cond
										 std::function<vecd (const vecd &, const vecd &, size_t, size_t)> wunc, // first argument is array of x_i, second — array of y(x_i,t); return value — array of w(x_i,y_i)
										 const vecd &x, // array with (non)uniform grid
										 vecd &y// array with initial condition and for results
									 	){
	const size_t N = std::min(x.size(), y.size()); // N_x+1
	auto W = wunc(x, y, 1, N-1);
	vecd K_0(N), K_1(N), frac(N), a(N), b(N), f(N);
	for ( size_t i = 1; i < N-1; ++i ){
		a[i] = 2. * ( x[i+1] - x[i] ) / ( x[i+1] - x[i-1] );
		b[i] = 2. * ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
		frac[i] = ( x[i+1] - x[i] ) * ( x[i] - x[i-1] ) / tau;
	}
    frac[N-1] = ( x[N-1] - x[N-2] ) * ( x[N-1] - x[N-2] ) * 0.5 / tau;
    for ( size_t i = 1; i < N; ++i ){
        f[i] = frac[i] * W[i];
		K_1[i] = frac[i] * W[i] / y[i];
    }

	vecd alpha(N), beta(N);
    double c;
    do {
		K_0 = K_1;
		alpha[1] = 0.;
		beta[1] = left_bounder_cond;
		for ( size_t i = 1; i < N-1; ++i ){
			c = 2. + K_1[i];
			alpha[i+1] = b[i] / ( c - alpha[i] * a[i] );
			beta[i+1] = ( beta[i] * a[i] + f[i] ) / ( c - alpha[i] * a[i] );
		}
		y[N-1] = ( (x[N-1] - x[N-2] ) * right_bounder_cond + f[N-1] + beta[N-1] ) / ( 1. + K_1[N-1] - alpha[N-1] );
		for ( size_t i = N-2; i > 0; --i ) {
			y[i] = alpha[i+1] * y[i+1] + beta[i+1];
		}
		y[0] = left_bounder_cond;
		W = wunc(x, y, 1, N-1);
		for ( size_t i = 1; i < N-1; ++i ){
			K_1[i] = frac[i] * W[i] / y[i];
		}
	} while( max_dif_rel(K_1, K_0, 1, N-2) > eps );
}
