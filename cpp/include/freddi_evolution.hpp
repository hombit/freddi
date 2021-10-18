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


#endif //FREDDI_FREDDI_EVOLUTION_HPP
