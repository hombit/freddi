#ifndef FREDDI_FREDDI_HPP
#define FREDDI_FREDDI_HPP

#include <functional>  // bind, function
#include <vector>

#include "arguments.hpp"
#include "spectrum.hpp"


class FreddiState;


class Freddi {
	typedef std::vector<double> vecd;
private:
	double Mdot_in_prev;
public:
	const double GM;
	const double eta;
	const double cosi;
	const double cosiOverD2;
	const OpacityRelated* oprel;
	std::function<vecd (const vecd&, const vecd&, unsigned int, unsigned int)> wunc;
	const FreddiArguments* args;
private:
	std::unique_ptr<FreddiState> state_;
private:
	FreddiState initializeState();
	void calculateRadialStructure();
	void truncateOuterRadius();
protected:
	vecd wunction(const vecd& h, const vecd& F, int first, int last) const;
	double Sigma_hot_disk(double r) const;
public:
	Freddi(const FreddiArguments& args);
	void step(double tau);
	inline void step() { return step(args->calc->tau); }
	std::vector<FreddiState> evolve();
public:
	inline const FreddiState& get_state() const { return *state_; }
};


class FreddiState {
	typedef std::vector<double> vecd;
private:
	const Freddi* freddi;
public:
	FreddiState(const Freddi* freddi);
	FreddiState(const FreddiState&) = default;
	FreddiState(FreddiState&&) = default;
	friend Freddi;
private:
	double Mdot_in = 0;
	double Mdot_out = 0;
	double Lx;
	double t = 0;
	size_t i_t = 0;
	unsigned int Nx;
	vecd h;
	vecd R;
	vecd F;
	vecd W;
	vecd Tph;
	vecd Tph_vis;
	vecd Tph_X;
	vecd Tirr;
	vecd Cirr;
	vecd Sigma;
	vecd Height;
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
	inline const vecd& get_Tph() const { return Tph; }
	inline const vecd& get_Tph_vis() const { return Tph_vis; }
	inline const vecd& get_Tirr() const { return Tirr; }
	inline const vecd& get_Cirr() const { return Cirr; }
	inline const vecd& get_Sigma() const { return Sigma; }
	inline const vecd& get_Height() const { return Height; }
	inline const double magnitude(double lambda, double F0) const {
		return -2.5 * log10( I_lambda(R, Tph, lambda) * freddi->cosiOverD2 / F0 );
	}
	inline const double flux(double lambda) const {
		return I_lambda(R, Tph, lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * freddi->cosiOverD2;
	}
	inline double mU() const { return magnitude(lambdaU, irr0U); }
	inline double mB() const { return magnitude(lambdaB, irr0B); }
	inline double mV() const { return magnitude(lambdaV, irr0V); }
	inline double mR() const { return magnitude(lambdaR, irr0R); }
	inline double mI() const { return magnitude(lambdaI, irr0I); }
	inline double mJ() const { return magnitude(lambdaJ, irr0J); }
	double integrate(const vecd& values) const;
	inline double Mdisk() const { return integrate(Sigma); }
private:
	void increase_t(double tau);
};

#endif //FREDDI_FREDDI_HPP
