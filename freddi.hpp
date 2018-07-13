#ifndef FREDDI_FREDDI_HPP
#define FREDDI_FREDDI_HPP

#include <functional>  // bind, function
#include <vector>

#include "arguments.hpp"
#include "spectrum.hpp"


class Freddi {
	typedef std::vector<double> vecd;
private:
	double Mdot_in = 0.;
	double Mdot_in_prev;
	double Mdot_out = 0.;
	double Lx = 0.;
	double t = 0.;
	unsigned int Nx;
	vecd h;
	vecd R;
	vecd F;
	vecd W;
	vecd Tph;
	vecd Tph_vis;
	vecd Tph_X;
	vecd Tirr;
	vecd Sigma;
	vecd Height;
public:
	const double GM;
	const double eta;
	const double h_in;
	const double h_out;
	const double cosi;
	const double cosiOverD2;
	const OpacityRelated* oprel;
	std::function<vecd (const vecd&, const vecd&, unsigned int, unsigned int)> wunc;
public:
	const FreddiArguments* args;
private:
	void initializeRadialStructure();
	void calculateRadialStructure();
	void truncateOuterRadius();
public:
	vecd wunction(const vecd& h, const vecd& F, int first, int last) const;
	double Sigma_hot_disk(double r) const;
public:
	Freddi(const FreddiArguments& args);
public:
	void next();
public:
	inline double get_Mdot_in() const { return Mdot_in; }
	inline double get_Mdot_out() const { return Mdot_out; }
	inline double get_Lx() const { return Lx; }
	inline double get_t() const { return t; }
	inline unsigned int get_Nx() const { return Nx; }
	inline const vecd& get_h() const { return h; }
	inline const vecd& get_R() const { return R; }
	inline const vecd& get_F() const { return F; }
	inline const vecd& get_W() const { return W; }
	inline const double magnitude(double lambda, double F0) const {
		return -2.5 * log10( I_lambda(R, Tph, lambda) * cosiOverD2 / F0 );
	}
	inline double mU() const { return magnitude(lambdaU, irr0U); }
	inline double mB() const { return magnitude(lambdaB, irr0B); }
	inline double mV() const { return magnitude(lambdaV, irr0V); }
	inline double mR() const { return magnitude(lambdaR, irr0R); }
	inline double mI() const { return magnitude(lambdaI, irr0I); }
	inline double mJ() const { return magnitude(lambdaJ, irr0J); }
	double integrate(const vecd& values) const;
	inline double Mdisk() const { return integrate(Sigma); }
};

void freddi(const FreddiArguments& args);

#endif //FREDDI_FREDDI_HPP
