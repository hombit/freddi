#include "nonlinear_diffusion.hpp"
#include "exceptions.hpp"  // gvlipunova added for DiscEqFailException, did not know where else....
#include <algorithm>  // std::max, std::min
#include <cmath>      // INFINITY
#include <cstdio>     // std::fprintf

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

void nonlinear_diffusion_nonuniform_wind_1_2_ (
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
	int iter_sol = 0;
	int maxiter = 100;
	int flag_F_negative = 0;
	do {
		iter_sol++;
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
			//if (i == last-1) {std::printf("ys: %e %e\n", y[i], y[i+1] );}
		}
		if (y[last]<0) {flag_F_negative = 1; }
		y[first] = left_bounder_cond;
		W = wunc(x, y, first + 1, last);
		for (size_t i = 1; i <= last - 1; ++i) {
			K_1[i] = frac[i] * W[i] / y[i];
		}
	} while ((max_dif_rel(K_1, K_0, 1, last - 1) > eps) && (iter_sol <=maxiter));
	 //std::printf ("flag_F_negative=%d  ", flag_F_negative);
	 if (iter_sol >= maxiter) { 
	     throw std::invalid_argument("Disc equation failed to converge. If you set --initialcond=gaussF, try move gausssigma and gaussmu parameters");
	     throw DiscEqFailException();
	} 
}

namespace {

// One implicit solve of the disc equation with the Picard iteration on K = W/y.
//
// For the iteration matrix of one Picard iteration the solution is affine in the outer
// boundary condition bc (= dF/dh at the outer edge): y_i = yA_i + yB_i * bc, with yA the
// solution for bc = 0 and yB the response to bc = 1. So the range of bc that keeps
// F >= 0 everywhere follows from one extra back-substitution, and when the requested
// bc lies outside it, bc is clamped to the bound (moved toward zero only, with a 0.1%
// margin so that no F is exactly zero). This replaces the former loop that lowered bc
// by 1% and re-solved, up to 100 times, without reporting the value it ended with.
//
// Returns true if the Picard iteration converged. n_clamped counts the iterations in
// which the clamp acted; bc_eff is the boundary condition of the final iteration.
// Throws DiscEqFailException if F would be negative even with bc = 0 (an interior
// problem that no boundary condition can fix).
bool solve_with_clamped_outer_bc(
        const double tau,
        const double eps,
        const double left_bounder_cond,
        const double bc_request,
        const vecd &A,
        const vecd &B,
        const vecd &C,
        const std::function<vecd (const vecd &, const vecd &, size_t, size_t)>& wunc,
        const vecd &x,
        vecd &y,
        size_t first, size_t last,
        int& n_clamped, double& bc_eff
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
        K_1[i] = f[i] / y[i];
    }
    K_1[last] = frac[last] * W[last] / y[last];

    vecd alpha(last + 1), beta(last + 1), yA(last + 1), yB(last + 1);
    const double dx_last = x[last] - x[last - 1];
    double c;
    int iter_sol = 0;
    const int maxiter = 100;
    n_clamped = 0;
    bc_eff = bc_request;

    do {
        iter_sol++;
        K_0 = K_1;
        alpha[first + 1] = 0.;
        beta[first + 1] = left_bounder_cond;
        for (size_t i = first + 1; i <= last - 1; ++i) {
            c = c0[i] + K_1[i];
            alpha[i + 1] = b[i] / (c - alpha[i] * a[i]);
            beta[i + 1] = (beta[i] * a[i] + f[i]) / (c - alpha[i] * a[i]);
        }

        // Affine dependence of the solution on the outer boundary condition:
        const double denom = c0[last] + K_1[last] - alpha[last] * a[last];
        yA[last] = (f[last] + beta[last] * a[last]) / denom;
        yB[last] = dx_last / denom;
        for (size_t i = last - 1; i > first; --i) {
            yA[i] = alpha[i + 1] * yA[i + 1] + beta[i + 1];
            yB[i] = alpha[i + 1] * yB[i + 1];
        }
        // Range of bc with y_i >= 0 for every i:
        double lower = -INFINITY, upper = INFINITY;
        for (size_t i = first + 1; i <= last; ++i) {
            if (yB[i] > 0.0) {
                lower = std::max(lower, -yA[i] / yB[i]);
            } else if (yB[i] < 0.0) {
                upper = std::min(upper, -yA[i] / yB[i]);
            } else if (yA[i] < 0.0) {
                lower = INFINITY;  // negative regardless of bc
            }
        }
        double bc = bc_request;
        if (lower > upper) {
            std::fprintf(stderr,
                "Disc equation: F would be negative for every outer boundary condition (interior problem, iteration %d)\n",
                iter_sol);
            throw DiscEqFailException();
        }
        if (bc < lower || bc > upper) {
            const double bound = (bc < lower) ? lower : upper;
            // Only move toward zero, as the former 1% loop did:
            if (bound != 0.0 && ((bound > 0.0) != (bc > 0.0))) {
                std::fprintf(stderr,
                    "Disc equation: F would be negative even with a zero outer boundary condition (requested %e, admissible from %e to %e)\n",
                    bc, lower, upper);
                throw DiscEqFailException();
            }
            bc = bound * (1.0 - 1e-3);
            n_clamped++;
        }
        bc_eff = bc;

        // Same back-substitution as before, so that unclamped solves give identical numbers:
        y[last] = (dx_last * bc + f[last] + beta[last] * a[last]) / denom;
        for (size_t i = last - 1; i > first; --i) {
            y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
        }
        y[first] = left_bounder_cond;
        W = wunc(x, y, first + 1, last);
        for (size_t i = 1; i <= last - 1; ++i) {
            K_1[i] = frac[i] * W[i] / y[i];
        }
    } while ((max_dif_rel(K_1, K_0, 1, last - 1) > eps) && (iter_sol <= maxiter));

    return iter_sol <= maxiter;
}

} // namespace

void nonlinear_diffusion_nonuniform_wind_1_2 (
        const double tau,
        const double eps,
        const double left_bounder_cond,
        const double right_bounder_cond,
        const vecd &A,
        const vecd &B,
        const vecd &C,
        const std::function<vecd (const vecd &, const vecd &, size_t, size_t)>& wunc,
        const vecd &x,
        vecd &y,
        size_t first, size_t last,
        int* bc_clamped, double* bc_factor
) {
    const vecd initial_y = y;  // restored before each attempt
    double bc_request = right_bounder_cond;
    const int max_attempts = 4;  // one solve plus three retries with the boundary condition halved

    for (int attempt = 1; ; ++attempt) {
        y = initial_y;
        int n_clamped = 0;
        double bc_eff = bc_request;
        const bool converged = solve_with_clamped_outer_bc(
            tau, eps, left_bounder_cond, bc_request, A, B, C, wunc, x, y, first, last,
            n_clamped, bc_eff);

        const double factor = (right_bounder_cond != 0.0) ? bc_eff / right_bounder_cond : 1.0;
        if (bc_clamped != nullptr) {
            *bc_clamped = std::max(*bc_clamped, n_clamped);  // per time step: worst solve
        }
        if (bc_factor != nullptr) {
            *bc_factor = std::min(*bc_factor, factor);
        }
        if (n_clamped > 0) {
            std::fprintf(stderr,
                "Disc equation: outer boundary condition clamped from %e to %e (factor %.4f) in %d iterations to keep F >= 0\n",
                right_bounder_cond, bc_eff, factor, n_clamped);
        }
        if (converged) {
            return;
        }
        if (attempt >= max_attempts) {
            std::fprintf(stderr,
                "Disc equation: nonlinear iteration did not converge in %d attempts (last outer boundary condition %e)\n",
                attempt, bc_request);
            throw DiscEqFailException();
        }
        bc_request *= 0.5;
        std::fprintf(stderr,
            "Disc equation: nonlinear iteration did not converge; retrying with the outer boundary condition halved to %e (attempt %d)\n",
            bc_request, attempt + 1);
    }
}
