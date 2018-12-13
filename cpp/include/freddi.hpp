#ifndef FREDDI_FREDDI_HPP
#define FREDDI_FREDDI_HPP

#include <functional>  // bind, function
#include <vector>

#include <boost/optional.hpp>

#include "arguments.hpp"
#include "spectrum.hpp"


class FreddiState;


class FreddiEvolution {
	typedef std::vector<double> vecd;
protected:
	double Mdot_in_prev;
public:
	const size_t Nt;
	const double GM;
	const double eta;
	const double cosi;
	const double cosiOverD2;
	const OpacityRelated* oprel;
	std::function<vecd (const vecd&, const vecd&, size_t, size_t)> wunc;
	const FreddiArguments* args;
protected:
	std::unique_ptr<FreddiState> state_;
	void truncateOuterRadius();
protected:
	vecd wunction(const vecd& h, const vecd& F, size_t first, size_t last) const;
	double Sigma_hot_disk(double r) const;
	inline void set_Mdot_in_prev(double Mdot_in) { Mdot_in_prev = Mdot_in; }
	void set_Mdot_in_prev();
public:
	FreddiEvolution(const FreddiArguments& args);
	virtual void step(double tau);
	inline void step() { return step(args->calc->tau); }
	std::vector<FreddiState> evolve();
	inline FreddiState& state() { return *state_; }
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
	const NeutronStarArguments* args_ns;
private:
	void truncateInnerRadius();
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
	using FreddiEvolution::step;
	virtual void step(double tau) override;
};


class FreddiState {
	typedef std::vector<double> vecd;
private:
	class BasicWind {
	protected:
		vecd A_, B_, C_;
		const FreddiState* state;
	public:
		BasicWind(const FreddiState* state);
		virtual void update(const FreddiState* _state) { state = _state; }
		inline const vecd& A() const { return A_; }
		inline const vecd& B() const { return B_; }
		inline const vecd& C() const { return C_; }
	};

	class NoWind: public BasicWind {
	public:
		NoWind(const FreddiState* state): BasicWind(state) {}
		~NoWind() = default;
	};

	class SS73CWind: public BasicWind {
	public:
		SS73CWind(const FreddiState* state);
		~SS73CWind() = default;
	};

	class Cambier2013Wind: public BasicWind {
	private:
		// windparams
		const double kC;
		const double R_IC2out;
	public:
		Cambier2013Wind(const FreddiState* state);
		~Cambier2013Wind() = default;
	};

	class __testA__Wind: public BasicWind {
	private:
		const double kA;
	public:
		__testA__Wind(const FreddiState* state);
		~__testA__Wind() = default;
	};

	class __testB__Wind: public BasicWind {
	private:
		const double kB;
	public:
		__testB__Wind(const FreddiState* state);
		~__testB__Wind() = default;
	};

	class __testC__Wind: public BasicWind {
	private:
		const double kC;
	public:
		__testC__Wind(const FreddiState* state);
		~__testC__Wind() = default;
	};

	class __testC_q0_Shields1986__: public BasicWind {
	private:
		const double kC;
		const double R_windmin2out;
	public:
		__testC_q0_Shields1986__(const FreddiState* state);
		~__testC_q0_Shields1986__() = default;
		virtual void update(const FreddiState* state) override;
	};
private:
	const FreddiEvolution* freddi;
protected:
	void initializeGrid();
	void initializeF();
	void initializeWind();
public:
	explicit FreddiState(const FreddiEvolution* freddi);
	explicit FreddiState(const FreddiNeutronStarEvolution* freddi);
	FreddiState(const FreddiState&, double tau);
	FreddiState(const FreddiState&) = default;
	FreddiState(FreddiState&&) = default;
	FreddiState& operator=(const FreddiState&) = default;
	FreddiState& operator=(FreddiState&&) = delete;
	friend FreddiEvolution;
	friend FreddiNeutronStarEvolution;
private:
	double Mdot_out_ = 0;
	double F_in_ = 0;
	double t_ = 0;
	size_t i_t_ = 0;
	size_t Nx_;
	vecd h_, R_;
	vecd F_;
public:
	inline double Mdot_in() const { return (F()[1] - F()[0]) / (h()[1] - h()[0]); }
	inline double Mdot_out() const { return Mdot_out_; }
	inline double F_in() const { return F_in_; }
	inline const vecd& h() const { return h_; }
	inline const vecd& R() const { return R_; }
	inline const vecd& F() const { return F_; }
	inline double t() const { return t_; }
	inline size_t i_t() const { return i_t_; };
	inline size_t Nx() const { return Nx_; }
private:
	boost::optional<double> Mdisk_;
	boost::optional<double> Lx_;
	boost::optional<double> mU_, mB_, mV_, mR_, mI_, mJ_;
	boost::optional<vecd> W_, Tph_, Qx_, Tph_vis_, Tph_X_, Tirr_, Cirr_, Sigma_, Height_;
	std::shared_ptr<BasicWind> wind_;
private:
	double lazy_magnitude(boost::optional<double>& m, double lambda, double F0);
	double lazy_integrate(boost::optional<double>& x, const vecd& values);
	const vecd& Tph_X();
	const vecd& Qx();
public:
	double Lx();
	const vecd& W();
	const vecd& Sigma();
	const vecd& Tph();
	const vecd& Tph_vis();
	const vecd& Tirr();
	const vecd& Cirr();
	const vecd& Height();
	inline double magnitude(double lambda, double F0) {
		return -2.5 * log10( I_lambda(R(), Tph(), lambda) * freddi->cosiOverD2 / F0 );
	}
	inline double flux(double lambda) {
		return I_lambda(R(), Tph(), lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * freddi->cosiOverD2;
	}
	inline double mU() { return lazy_magnitude(mU_, lambdaU, irr0U); }
	inline double mB() { return lazy_magnitude(mB_, lambdaB, irr0B); }
	inline double mV() { return lazy_magnitude(mV_, lambdaV, irr0V); }
	inline double mR() { return lazy_magnitude(mR_, lambdaR, irr0R); }
	inline double mI() { return lazy_magnitude(mI_, lambdaI, irr0I); }
	inline double mJ() { return lazy_magnitude(mJ_, lambdaJ, irr0J); }
	double integrate(const vecd& values) const;
	inline double Mdisk() { return lazy_integrate(Mdisk_, Sigma()); }
};

#endif //FREDDI_FREDDI_HPP
