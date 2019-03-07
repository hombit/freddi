#ifndef FREDDI_FREDDI_EVOLUTION_HPP
#define FREDDI_FREDDI_EVOLUTION_HPP

#include <functional>  // bind, function
#include <vector>

#include <boost/optional.hpp>

#include "arguments.hpp"
#include "freddi_state.hpp"
#include "spectrum.hpp"


class FreddiEvolution: public FreddiState {
protected:
	virtual void truncateOuterRadius();
	virtual void truncateInnerRadius() {}
protected:
	virtual vecd wunction(const vecd& h, const vecd& F, size_t first, size_t last) const;
public:
	FreddiEvolution(const FreddiArguments& args);
	FreddiEvolution(const FreddiEvolution&) = delete;
	virtual void step(double tau);
	inline void step() { return step(args().calc->tau); }
};


class FreddiNeutronStarEvolution: public FreddiEvolution {
public:
	const NeutronStarArguments* args_ns;
	const double k_t = 1. / 3.;
	const double xi = 0.7;
	const double xi_pow_minus_7_2;
	const double R_m_min;
	const double mu_magn;
	const double R_dead;
	const double F_dead;
	const double R_cor;
	const double inverse_beta;
public:
	const vecd Fmagn;
	const vecd dFmagn_dh;
	const vecd d2Fmagn_dh2;
protected:
	const vecd initialize_Fmagn() const;
	const vecd initialize_dFmagn_dh() const;
	const vecd initialize_d2Fmagn_dh2() const;
public:
public:
	using FreddiEvolution::step;
protected:
	virtual void truncateInnerRadius() override;
	virtual double Mdot_in() const override;
	virtual vecd windC() const override;
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
	FreddiNeutronStarEvolution(const FreddiNeutronStarEvolution&) = delete;
};


#endif //FREDDI_FREDDI_EVOLUTION_HPP
