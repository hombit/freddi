#ifndef FREDDI_FREDDI_EVOLUTION_HPP
#define FREDDI_FREDDI_EVOLUTION_HPP

#include <functional>  // bind, function
#include <iterator>
#include <vector>

#include <boost/optional.hpp>

#include "arguments.hpp"
#include "freddi_state.hpp"
#include "spectrum.hpp"


template<typename T>
class EvolutionIterator: public std::iterator<std::forward_iterator_tag, size_t, ptrdiff_t, const T*, T&> {
private:
	T* evolution;
	size_t i_t;
public:
	EvolutionIterator(T* evolution): evolution(evolution), i_t(evolution->i_t()) {}
	EvolutionIterator(size_t i_t): evolution(nullptr), i_t(i_t) {}
	EvolutionIterator& operator++() {
		i_t++;
		return *this;
	}
	EvolutionIterator operator++(int) {
		auto tmp = *this;
		++(*this);
		return tmp;
	}
	bool operator==(EvolutionIterator other) const { return i_t == other.i_t; }
	bool operator!=(EvolutionIterator other) const { return i_t != other.i_t; }
	T& operator*() const {
		if (evolution->i_t() < i_t) {
			evolution->step();
		}
		return *evolution;
	}
};


class FreddiEvolution: public FreddiState {
protected:
	virtual void truncateOuterRadius();
	virtual void truncateInnerRadius() {}
protected:
	virtual vecd wunction(const vecd& h, const vecd& F, size_t first, size_t last) const;
public:
	FreddiEvolution(const FreddiArguments& args);
	explicit FreddiEvolution(const FreddiEvolution&) = default;
	virtual void step(double tau);
	inline void step() { return step(args().calc->tau); }
public:
	using iterator = EvolutionIterator<FreddiEvolution>;
	inline iterator begin() { return {this}; }
	inline iterator end() { return {Nt() + 1}; }
};


class FreddiNeutronStarEvolution: public FreddiEvolution {
private:
	class MagneticFieldStructure {
	public:
		const NeutronStarArguments args_ns;
		const double k_t = 1. / 3.;
		const double xi = 0.7;
		const double xi_pow_minus_7_2;
		const double R_m_min;
		const double mu_magn;
		const double R_dead;
		const double F_dead;
		const double R_cor;
		const double inverse_beta;
		const double epsilon_Alfven;
		const vecd Fmagn;
		const vecd dFmagn_dh;
		const vecd d2Fmagn_dh2;
		MagneticFieldStructure(const NeutronStarArguments& args_ns, FreddiEvolution* evolution);
	protected:
		vecd initialize_Fmagn(FreddiEvolution* evolution) const;
		vecd initialize_dFmagn_dh(FreddiEvolution* evolution) const;
		vecd initialize_d2Fmagn_dh2(FreddiEvolution* evolution) const;
	};
private:
	std::shared_ptr<MagneticFieldStructure> magn_field_str_;
public:
	inline double k_t() const { return magn_field_str_->k_t; }
	inline double xi() const { return magn_field_str_->xi; }
	inline double xi_pow_minus_7_2() const { return magn_field_str_->xi_pow_minus_7_2; }
	inline double R_m_min() const { return magn_field_str_->R_m_min; }
	inline double mu_magn() const { return magn_field_str_->mu_magn; }
	inline double R_dead() const { return magn_field_str_->R_dead; }
	inline double F_dead() const { return magn_field_str_->F_dead; }
	inline double R_cor() const { return magn_field_str_->R_cor; }
	inline double inverse_beta() const { return magn_field_str_->inverse_beta; }
	inline double epsilon_Alfven() const { return magn_field_str_->epsilon_Alfven; }
	inline const vecd& Fmagn() const { return magn_field_str_->Fmagn; }
	inline const vecd& dFmagn_dh() const { return magn_field_str_->dFmagn_dh; }
	inline const vecd& d2Fmagn_dh2() const { return magn_field_str_->d2Fmagn_dh2; }
public:
	using FreddiEvolution::step;
protected:
	virtual void truncateInnerRadius() override;
	virtual double Mdot_in() const override;
	virtual vecd windC() const override;
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
	explicit FreddiNeutronStarEvolution(const FreddiNeutronStarEvolution&) = default;
public:
	using iterator = EvolutionIterator<FreddiNeutronStarEvolution>;
	inline iterator begin() { return {this}; }
	inline iterator end() { return {Nt() + 1}; }
};


#endif //FREDDI_FREDDI_EVOLUTION_HPP
