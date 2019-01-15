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
	virtual void step(double tau);
	inline void step() { return step(args.calc->tau); }
	inline const FreddiState& state() { return static_cast<FreddiState&>(*this); }
};


class FreddiNeutronStarEvolution: public FreddiEvolution {
public:
	const double k_t = 1. / 3.;
	const double xi = 0.7;
	const double R_m_min;
	const double mu_magn;
	const double R_dead;
	const double R_cor;
	const double xi_pow_minus_7_2;
	const double inverse_beta;
	const NeutronStarArguments* args_ns;
private:
	vecd magnetic_windC;
public:
	using FreddiEvolution::step;
protected:
	virtual void truncateInnerRadius() override;
	virtual const vecd windC() const override;
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
};


#endif //FREDDI_FREDDI_EVOLUTION_HPP
