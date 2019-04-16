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
protected:
private:
	class NeutronStarStructure {
	public:
		NeutronStarArguments args_ns;
		double k_t = 1. / 3.;
		double xi = 0.7;
		double xi_pow_minus_7_2;
		double R_x;
		double R_m_min;
		double mu_magn;
		double R_dead;
		double F_dead;
		double R_cor;
		double inverse_beta;
		double epsilon_Alfven;
		double hot_spot_area;
		vecd Fmagn;
		vecd dFmagn_dh;
		vecd d2Fmagn_dh2;
		NeutronStarStructure(const NeutronStarArguments& args_ns, FreddiEvolution* evolution);
	protected:
		vecd initialize_Fmagn(FreddiEvolution* evolution) const;
		vecd initialize_dFmagn_dh(FreddiEvolution* evolution) const;
		vecd initialize_d2Fmagn_dh2(FreddiEvolution* evolution) const;
	};

	struct NeutronStarOptionalStructure {
		boost::optional<double> Lx_ns;
	};

	class BasicNSMdotFraction {
	public:
		BasicNSMdotFraction() = default;
		virtual ~BasicNSMdotFraction() = 0;
		virtual double operator()(double R_to_Rcor) = 0;
	};

	class NoOutflowNSMdotFraction: public BasicNSMdotFraction {
	public:
		NoOutflowNSMdotFraction() = default;
		~NoOutflowNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};

	class PropellerNSMdotFraction: public BasicNSMdotFraction {
	public:
		PropellerNSMdotFraction() = default;
		~PropellerNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};

	class CorotationBlockNSMdotFraction: public BasicNSMdotFraction {
	public:
		CorotationBlockNSMdotFraction() = default;
		~CorotationBlockNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};

	// https://arxiv.org/pdf/1010.1528.pdf Eksi-Kutlu (2010)
	class EksiKultu2010NSMdotFraction: public BasicNSMdotFraction {
	public:
		EksiKultu2010NSMdotFraction() = default;
		~EksiKultu2010NSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};

	class Romanova2018NSMdotFraction: public BasicNSMdotFraction {
	private:
		const double par1;
		const double par2;
	public:
		Romanova2018NSMdotFraction(double par1, double par2);
		~Romanova2018NSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};

	class GeometricalNSMdotFraction: public BasicNSMdotFraction {
	private:
		const double chi;
	public:
		GeometricalNSMdotFraction(double chi);
		~GeometricalNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) override;
	};
private:
	std::shared_ptr<const NeutronStarStructure> ns_str_;
	NeutronStarOptionalStructure ns_opt_str_;
	std::shared_ptr<BasicNSMdotFraction> fp_;
private:
	void initializeNsMdotFraction();
// ns_str_
public:
	inline double k_t() const { return ns_str_->k_t; }
	inline double xi() const { return ns_str_->xi; }
	inline double xi_pow_minus_7_2() const { return ns_str_->xi_pow_minus_7_2; }
	inline double R_x() const { return ns_str_->R_x; }
	inline double R_m_min() const { return ns_str_->R_m_min; }
	inline double mu_magn() const { return ns_str_->mu_magn; }
	inline double R_dead() const { return ns_str_->R_dead; }
	inline double F_dead() const { return ns_str_->F_dead; }
	inline double R_cor() const { return ns_str_->R_cor; }
	inline double inverse_beta() const { return ns_str_->inverse_beta; }
	inline double epsilon_Alfven() const { return ns_str_->epsilon_Alfven; }
	inline double hot_spot_area() const { return ns_str_->hot_spot_area; }
	inline const vecd& Fmagn() const { return ns_str_->Fmagn; }
	inline const vecd& dFmagn_dh() const { return ns_str_->dFmagn_dh; }
	inline const vecd& d2Fmagn_dh2() const { return ns_str_->d2Fmagn_dh2; }
// ns_opt_str_
public:
	double Lbol_ns() const;
	double T_hot_spot() const;
	double Lx_ns();
// fp_
public:
	inline double fp(double R) const { return (*fp_)(R / R_cor()); }
	inline double fp() const { return (*fp_)(R()[first()]); }
public:
	using FreddiEvolution::step;
protected:
	virtual void invalidate_optional_structure() override;
	virtual void truncateInnerRadius() override;
	virtual double Mdot_in() const override;
	virtual vecd windC() const override;
	virtual const vecd& Qx() override;
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
	explicit FreddiNeutronStarEvolution(const FreddiNeutronStarEvolution&) = default;
public:
	using iterator = EvolutionIterator<FreddiNeutronStarEvolution>;
	inline iterator begin() { return {this}; }
	inline iterator end() { return {Nt() + 1}; }
public:
	virtual double eta_ns() const;
};


#endif //FREDDI_FREDDI_EVOLUTION_HPP
