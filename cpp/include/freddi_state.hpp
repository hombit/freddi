#ifndef FREDDI_FREDDI_STATE_HPP
#define FREDDI_FREDDI_STATE_HPP

#include <functional>  // bind, function
#include <vector>

#include <boost/optional.hpp>

#include "arguments.hpp"
#include "spectrum.hpp"
#include "util.hpp"


class FreddiState {
protected:
	typedef std::function<vecd (const vecd&, const vecd&, size_t, size_t)> wunc_t;
private:
	class BasicWind {
	protected:
		vecd A_, B_, C_;
	public:
		explicit BasicWind(const FreddiState&);
		virtual ~BasicWind() = 0;
		virtual BasicWind* clone() const = 0;
		virtual void update(const FreddiState&) {}
		inline const vecd& A() const { return A_; }
		inline const vecd& B() const { return B_; }
		inline const vecd& C() const { return C_; }
	};

	class NoWind: public BasicWind {
	public:
		explicit NoWind(const FreddiState& state): BasicWind(state) {}
		~NoWind() override = default;
		NoWind(const NoWind&) = default;
		virtual NoWind* clone() const override { return new NoWind(*this); }
	};

	class SS73CWind: public BasicWind {
	public:
		explicit SS73CWind(const FreddiState& state);
		~SS73CWind() override = default;
		SS73CWind(const SS73CWind&) = default;
		virtual SS73CWind* clone() const override { return new SS73CWind(*this); }
	};

	class Cambier2013Wind: public BasicWind {
	private:
		// windparams
		const double kC;
		const double R_IC2out;
	public:
		explicit Cambier2013Wind(const FreddiState& state);
		~Cambier2013Wind() override = default;
		Cambier2013Wind(const Cambier2013Wind&) = default;
		virtual Cambier2013Wind* clone() const override { return new Cambier2013Wind(*this); }
	};

	class __testA__Wind: public BasicWind {
	private:
		// windparams
		const double kA;
	public:
		explicit __testA__Wind(const FreddiState& state);
		~__testA__Wind() override = default;
		__testA__Wind(const __testA__Wind&) = default;
		virtual __testA__Wind* clone() const override { return new __testA__Wind(*this); }
	};

	class __testB__Wind: public BasicWind {
	private:
		// windparams
		const double kB;
	public:
		explicit __testB__Wind(const FreddiState& state);
		~__testB__Wind() override = default;
		__testB__Wind(const __testB__Wind&) = default;
		virtual __testB__Wind* clone() const override { return new __testB__Wind(*this); }
	};

	class __testC__Wind: public BasicWind {
	private:
		// windparams
		const double kC;
	public:
		explicit __testC__Wind(const FreddiState& state);
		~__testC__Wind() override = default;
		__testC__Wind(const __testC__Wind&) = default;
		virtual __testC__Wind* clone() const override { return new __testC__Wind(*this); }
	};

	class __testC_q0_Shields1986__: public BasicWind {
	private:
		// windparams
		const double kC;
		const double R_windmin2out;
	public:
		explicit __testC_q0_Shields1986__(const FreddiState& state);
		~__testC_q0_Shields1986__() override = default;
		__testC_q0_Shields1986__(const __testC_q0_Shields1986__&) = default;
		virtual __testC_q0_Shields1986__* clone() const override { return new __testC_q0_Shields1986__(*this); }
		virtual void update(const FreddiState&) override;
	};

	class __Unstedy_Test_Hunter__: public BasicWind {
	private:
		// windparams
		const double fXI;
		const double T_iC;
	public:
		explicit __Unstedy_Test_Hunter__(const FreddiState& state);
		~__Unstedy_Test_Hunter__() override = default;
		__Unstedy_Test_Hunter__(const __Unstedy_Test_Hunter__&) = default;
		virtual __Unstedy_Test_Hunter__* clone() const override { return new __Unstedy_Test_Hunter__(*this); }
		virtual void update(const FreddiState&) override;
	};

	class DiskStructure {
	public:
		FreddiArguments args;
		size_t Nt;
		size_t Nx;
		double GM;
		double eta;
		double cosi;
		double distance;
		double cosiOverD2;
		OpacityRelated oprel;
		vecd h;
		vecd R;
		wunc_t wunc;
		DiskStructure(const FreddiArguments& args, const wunc_t& wunc);
	private:
		static vecd initialize_h(const FreddiArguments& args, size_t Nx);
		static vecd initialize_R(const vecd& h, double GM);
	};

	class CurrentState {
	public:
		double F_in = 0;
		double Mdot_out;
		double Mdot_in_prev = -INFINITY;
		double t = 0;
		size_t i_t = 0;
		size_t first = 0;
		size_t last;
		vecd F;
		explicit CurrentState(const DiskStructure& str);
		CurrentState(const CurrentState&) = default;
		CurrentState& operator=(const CurrentState&) = default;
	private:
		static vecd initializeF(const DiskStructure& str);
	};

	struct DiskOptionalStructure {
		boost::optional<double> Mdisk;
		boost::optional<double> Lx;
		boost::optional<double> mU, mB, mV, mR, mI, mJ;
		boost::optional<vecd> W, Tph, Qx, Tph_vis, Tph_X, Tirr, Cirr, Sigma, Height;
	};

protected:
	std::shared_ptr<const DiskStructure> str_;
	CurrentState current_;
	DiskOptionalStructure opt_str_;
	std::unique_ptr<BasicWind> wind_;
public:
	FreddiState(const FreddiArguments& args, const wunc_t& wunc);
	FreddiState(const FreddiState&);
	FreddiState(FreddiState&&) = delete;
	FreddiState& operator=(const FreddiState&) = delete;
	FreddiState& operator=(FreddiState&&) = delete;
	virtual void step(double tau);
private:
	void initializeWind();
// str_
public:
	inline size_t Nt() const { return str_->Nt; }
	inline size_t Nx() const { return str_->Nx; }
	inline double GM() const { return str_->GM; }
	inline double eta() const { return str_->eta; }
	inline double cosi() const { return str_->cosi; }
	inline double distance() const { return str_->distance; }
	inline double cosiOverD2() const { return str_->cosiOverD2; }
	inline const OpacityRelated& oprel() const { return str_->oprel; }
	inline const wunc_t& wunc() const { return str_->wunc; }
	inline const FreddiArguments& args() const { return str_->args; }
	inline const vecd& h() const { return str_->h; }
	inline const vecd& R() const { return str_->R; }
// current_
public:
	inline double Mdot_out() const { return current_.Mdot_out; }
	inline double F_in() const { return current_.F_in; }
	inline const vecd& F() const { return current_.F; }
	inline double t() const { return current_.t; }
	inline size_t i_t() const { return current_.i_t; };
	inline size_t first() const { return current_.first; }
	inline size_t last() const { return current_.last; }
	inline double Mdot_in_prev() const { return current_.Mdot_in_prev; }
	inline void set_Mdot_in_prev(double Mdot_in) { current_.Mdot_in_prev = Mdot_in; }
	inline void set_Mdot_in_prev() { set_Mdot_in_prev(Mdot_in()); }
public:
	virtual double Mdot_in() const;
// wind_
protected:
	virtual vecd windA() const { return wind_->A(); }
	virtual vecd windB() const { return wind_->B(); }
	virtual vecd windC() const { return wind_->C(); }
// opt_str_
protected:
	inline void invalidate_disk_optional_structure() { opt_str_ = DiskOptionalStructure(); };
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
	double I_lambda(double lambda);
	double Luminosity(const vecd& T, double nu1, double nu2) const;
	inline double magnitude(const double lambda, const double F0) {
		return -2.5 * std::log10(I_lambda(lambda) * cosiOverD2() / F0);
	}
	inline double flux(const double lambda) {
		return I_lambda(lambda) * lambda*lambda / GSL_CONST_CGSM_SPEED_OF_LIGHT * cosiOverD2();
	}
	inline double mU() { return lazy_magnitude(opt_str_.mU, lambdaU, irr0U); }
	inline double mB() { return lazy_magnitude(opt_str_.mB, lambdaB, irr0B); }
	inline double mV() { return lazy_magnitude(opt_str_.mV, lambdaV, irr0V); }
	inline double mR() { return lazy_magnitude(opt_str_.mR, lambdaR, irr0R); }
	inline double mI() { return lazy_magnitude(opt_str_.mI, lambdaI, irr0I); }
	inline double mJ() { return lazy_magnitude(opt_str_.mJ, lambdaJ, irr0J); }
	inline double integrate(const vecd& values) const { return disk_radial_trapz(R(), values, first(), last()); }
	inline double integrate(std::function<double (size_t)> f) const { return disk_radial_trapz(R(), f, first(), last()); }
	inline double Mdisk() { return lazy_integrate(opt_str_.Mdisk, Sigma()); }
};

#endif //FREDDI_FREDDI_STATE_HPP
