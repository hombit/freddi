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



// \frac{dw}{dt}=\frac{d^2y}{dx^2} + A\frac{dy}{dx} + By + C, y=y(x,t) — ?, w = w (x,y)

void nonlinear_diffusion_nonuniform_wind_1_2 (
		const double tau,
		const double eps, // relative error for w
		const double left_bounder_cond, // y(left_border,Time+tau) = left_bounder_cond
		const double right_bounder_cond, // \frac{y(right_border,Time+tau)}{dx} = right_bounder_cond
		const vecd &A,
		const vecd &B,
		const vecd &C,
		const std::function<vecd (const vecd &, const vecd &, size_t, size_t)>& wunc, // first argument is array of x_i, second — array of y(x_i,t); return value — array of w(x_i,y_i)
		const vecd &x, // array with (non)uniform grid
		vecd &y, // array with initial condition and for results
		size_t first, size_t last // indexes of front and back elements
) {
	auto W = wunc(x, y, first + 1, last);
	vecd K_0(last + 1), K_1(last + 1), frac(last + 1), a(last + 1), b(last + 1), c0(last + 1), f(last + 1);
	for (size_t i = first + 1; i <= last - 1; ++i) {
		a[i] = (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]) * (2.0 - A[i] * (x[i + 1] - x[i]));
		b[i] = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]) * (2.0 + A[i] * (x[i] - x[i - 1]));
		c0[i] = 2.0 - A[i] * (x[i + 1] - 2 * x[i] + x[i - 1]) - B[i] * (x[i + 1] - x[i]) * (x[i] - x[i - 1]);
		frac[i] = (x[i + 1] - x[i]) * (x[i] - x[i - 1]) / tau;
	}
	a[last] = 1 - 0.5 * A[last] * (x[last] - x[last - 1]);
	c0[last] = a[last] - 0.5 * B[last] * (x[last] - x[last - 1]) * (x[last] - x[last - 1]);
	frac[last] = (x[last] - x[last - 1]) * (x[last] - x[last - 1]) * 0.5 / tau;
	for (size_t i = first + 1; i <= last; ++i) {
		f[i] = frac[i] * (W[i] + tau * C[i]);
	}
	for (size_t i = first + 1; i <= last - 1; ++i) {
//		K_1[i] = (f[i] + a[i] * y[i - 1] - c0[i] * y[i] + b[i] * y[i + 1]) / y[i];
//		K_1[i] = frac[i] * W[i] / y[i];
		K_1[i] = f[i] / y[i];
	}
//	K_1[last] = (f[last] + a[last] * y[last - 1] - c0[last] * y[last] + right_bounder_cond * (x[last] - x[last - 1])) / y[last];
//	K_1[last] = (f[last]) / y[last];
	K_1[last] = frac[last] * W[last] / y[last];

	vecd alpha(last + 1), beta(last + 1);
	double c;
	do {
		K_0 = K_1;
		alpha[first + 1] = 0.;
		beta[first + 1] = left_bounder_cond;
		for (size_t i = first + 1; i <= last - 1; ++i) {
			c = c0[i] + K_1[i];
			alpha[i + 1] = b[i] / (c - alpha[i] * a[i]);
			beta[i + 1] = (beta[i] * a[i] + f[i]) / (c - alpha[i] * a[i]);
		}
		y[last] = ((x[last] - x[last - 1]) * right_bounder_cond + f[last] + beta[last] * a[last]) /
				   (c0[last] + K_1[last] - alpha[last] * a[last]);
		for (size_t i = last - 1; i > first; --i) {
			y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
		}
		y[first] = left_bounder_cond;
		W = wunc(x, y, first + 1, last);
		for (size_t i = 1; i <= last - 1; ++i) {
			K_1[i] = frac[i] * W[i] / y[i];
		}
	} while (max_dif_rel(K_1, K_0, 1, last - 1) > eps);
}
