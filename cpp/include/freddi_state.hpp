#ifndef FREDDI_FREDDI_STATE_HPP
#define FREDDI_FREDDI_STATE_HPP

#include <functional>  // bind, function
#include <vector>

#include <boost/optional.hpp>

#include <arguments.hpp>
#include <passband.hpp>
#include <rochelobe.hpp>
#include <spectrum.hpp>
#include <star.hpp>
#include <util.hpp>


class FreddiState {
protected:
	typedef std::function<vecd (const vecd&, const vecd&, size_t, size_t)> wunc_t;
public:
	enum DiskIntegrationRegion {
		HotRegion,
		ColdRegion,
	};
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

	class testAWind: public BasicWind {
	private:
		// windparams
		const double kA;
	public:
		explicit testAWind(const FreddiState& state);
		~testAWind() override = default;
		testAWind(const testAWind&) = default;
		virtual testAWind* clone() const override { return new testAWind(*this); }
	};

	class testBWind: public BasicWind {
	private:
		// windparams
		const double kB;
	public:
		explicit testBWind(const FreddiState& state);
		~testBWind() override = default;
		testBWind(const testBWind&) = default;
		virtual testBWind* clone() const override { return new testBWind(*this); }
	};

	class testCWind: public BasicWind {
	private:
		// windparams
		const double kC;
	public:
		explicit testCWind(const FreddiState& state);
		~testCWind() override = default;
		testCWind(const testCWind&) = default;
		virtual testCWind* clone() const override { return new testCWind(*this); }
	};

	class testCq0Shields1986Wind: public BasicWind {
	private:
		// windparams
		const double kC;
		const double R_windmin2out;
	public:
		explicit testCq0Shields1986Wind(const FreddiState& state);
		~testCq0Shields1986Wind() override = default;
		testCq0Shields1986Wind(const testCq0Shields1986Wind&) = default;
		virtual testCq0Shields1986Wind* clone() const override { return new testCq0Shields1986Wind(*this); }
		virtual void update(const FreddiState&) override;
	};

	class Shields1986Wind: public BasicWind {
	//G. A. Shields, C. F. McKee, D. N. C. Lin, and M. C. Begelman. Compton-heated winds andcoronae above accretion disks. II - Instability and oscillations. ApJ, 306:90–106, July 1986.
	//doi:10.1086/164322
	private:
		// windparams
		const double f_X;
		const double X_f;
		const double T_iC;
	public:
		explicit Shields1986Wind(const FreddiState& state);
		~Shields1986Wind() override = default;
		Shields1986Wind(const Shields1986Wind&) = default;
		virtual Shields1986Wind* clone() const override { return new Shields1986Wind(*this); }
		virtual void update(const FreddiState&) override;
	};

	class Janiuk2015Wind: public BasicWind {
	//Janiuk A., Grzedzielski M., Capitanio F., Bianchi S., 2015, Interplay between heartbeat oscillations and wind outflow in microquasar IGR J17091-3624 A&A, 574, A92
	//doi:10.1051/0004-6361/201425003
	private:
		// windparams
		const double A_0;
		const double B_1;
	public:
		explicit Janiuk2015Wind(const FreddiState& state);
		~Janiuk2015Wind() override = default;
		Janiuk2015Wind(const  Janiuk2015Wind&) = default;
		virtual Janiuk2015Wind* clone() const override { return new Janiuk2015Wind(*this); }
		virtual void update(const FreddiState&) override;
	};

	class Woods1996Wind: public BasicWind {
	//D. T. Woods, R. I. Klein, J. I. Castor, C. F. McKee, and J. B. Bell. X-Ray–heated Coronae andWinds from Accretion Disks: Time-dependent Two-dimensional Hydrodynamics with AdaptiveMesh Refinement. ApJ, 461:767, April 1996
	//doi:10.1086/177101
	private:
		// windparams
		const double C_0;
		const double T_iC;
	public:
		explicit Woods1996Wind(const FreddiState& state);
		~Woods1996Wind() override = default;
		Woods1996Wind(const Woods1996Wind&) = default;
		virtual Woods1996Wind* clone() const override { return new Woods1996Wind(*this); }
		virtual void update(const FreddiState&) override;
	};

	class Woods1996ShieldsApproxWind : public BasicWind {
	//D. T. Woods, R. I. Klein, J. I. Castor, C. F. McKee, and J. B. Bell. X-Ray–heated Coronae andWinds from Accretion Disks: Time-dependent Two-dimensional Hydrodynamics with AdaptiveMesh Refinement. ApJ, 461:767, April 1996
	//doi:10.1086/177101
	private:
		// windparams
		const double Xi_max;
		const double T_iC;
		const double W_pow;
	public:
		explicit Woods1996ShieldsApproxWind(const FreddiState& state);
		~Woods1996ShieldsApproxWind() override = default;
		Woods1996ShieldsApproxWind(const Woods1996ShieldsApproxWind&) = default;
		virtual Woods1996ShieldsApproxWind* clone() const override { return new Woods1996ShieldsApproxWind(*this); }
		virtual void update(const FreddiState&) override;
	};
	
	class PeriodPaperWind : public BasicWind {
	// Avakyan, 2021
	// https://ui.adsabs.harvard.edu/abs/2021arXiv210511974A/
	private:
		// windparams
		const double W_pow;
	public:
		explicit PeriodPaperWind(const FreddiState& state);
		~PeriodPaperWind() override = default;
		PeriodPaperWind(const PeriodPaperWind&) = default;
		virtual PeriodPaperWind* clone() const override { return new PeriodPaperWind(*this); }
		virtual void update(const FreddiState&) override;
	};
	
protected:
	class BasicFreddiIrradiationSource {
	public:
		virtual ~BasicFreddiIrradiationSource() = 0;
		virtual double angular_dist(double mu) const = 0; // mu = cos(angle between ray and normal)
		virtual std::unique_ptr<IrrSource> irr_source(FreddiState& state, double luminosity) const = 0;
	protected:
		Vec3 position(const FreddiState& state) const;
		double Height2R(FreddiState& state) const;
	};

	class IsotropicFreddiIrradiationSource: public BasicFreddiIrradiationSource {
	public:
		~IsotropicFreddiIrradiationSource() override = default;
		double angular_dist(double mu) const override;
		std::unique_ptr<IrrSource> irr_source(FreddiState& state, double luminosity) const override;
	};

	class PlaneFreddiIrradiationSource: public BasicFreddiIrradiationSource {
	protected:
		const UnitVec3 normal;
	public:
		PlaneFreddiIrradiationSource();
		~PlaneFreddiIrradiationSource() override = default;
		double angular_dist(double mu) const override;
		std::unique_ptr<IrrSource> irr_source(FreddiState& state, double luminosity) const override;
	};

private:
	class DiskStructure {
	public:
		FreddiArguments args;
		size_t Nt;
		size_t Nx;
		double GM;
		double R_g;
		double eta;
		double semiaxis;
		double inclination;
		double cosi;
		double distance;
		double cosiOverD2;
		OpacityRelated oprel;
		vecd h;
		vecd R;
		wunc_t wunc;
	private:
		static vecd initialize_h(const FreddiArguments& args, size_t Nx);
		static vecd initialize_R(const vecd& h, double GM);
	public:
		DiskStructure(const FreddiArguments& args, const wunc_t& wunc);
	};

	class CurrentState {
	public:
		double Mdot_out;
		double Mdot_in_prev = -INFINITY;
		double t;
		size_t i_t;
		size_t first;
		size_t last;
		vecd F;
		double F_in;
		explicit CurrentState(const DiskStructure& str);
		CurrentState(const CurrentState&) = default;
		CurrentState& operator=(const CurrentState&) = default;
	private:
		static vecd initializeF(const DiskStructure& str);
	};

	struct DiskOptionalStructure {
		boost::optional<double> Mdisk;
		boost::optional<double> Lx;
		boost::optional<double> Mdot_wind;
		boost::optional<double> mU, mB, mV, mR, mI, mJ;
		boost::optional<vecd> W, Tph, Qx, Tph_vis, Tph_X, Tirr, Kirr, Sigma, Height;
	};

protected:
	std::shared_ptr<const DiskStructure> str_;
	CurrentState current_;
	DiskOptionalStructure opt_str_;
	std::unique_ptr<BasicWind> wind_;
	std::shared_ptr<BasicFreddiIrradiationSource> disk_irr_source_;
	RocheLobe star_roche_lobe_;
	IrradiatedStar star_;
public:
	FreddiState(const FreddiArguments& args, const wunc_t& wunc);
	explicit FreddiState(const FreddiState&);
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
	inline double R_g() const { return str_->R_g; }
	inline double eta() const { return str_->eta; }
	inline double semiaxis() const { return str_->semiaxis; }
	inline double inclination() const { return str_->inclination; }
	inline double cosi() const { return str_->cosi; }
	inline double distance() const { return str_->distance; }
	inline double cosiOverD2() const { return str_->cosiOverD2; }
	inline const OpacityRelated& oprel() const { return str_->oprel; }
	inline const wunc_t& wunc() const { return str_->wunc; }
	inline const FreddiArguments& args() const { return str_->args; }
	inline const vecd& h() const { return str_->h; }
	inline const vecd& R() const { return str_->R; }
	inline const vecd& lambdas() const { return str_->args.flux->lambdas; }
	inline Star& star() { return star_; }
	void replaceArgs(const FreddiArguments& args);  // Danger!
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
protected:
	inline void set_Mdot_in_prev(double Mdot_in) { current_.Mdot_in_prev = Mdot_in; }
	inline void set_Mdot_in_prev() { set_Mdot_in_prev(Mdot_in()); }
	virtual IrradiatedStar::sources_t star_irr_sources();
public:
	inline double omega_R(double r) const { return std::sqrt(GM() / (r*r*r)); }
	inline double omega_i(size_t i) const { return omega_R(R()[i]); }
	virtual double Mdot_in() const;
	virtual double Lbol_disk() const;
	double phase_opt() const;
// wind_
public:
	virtual vecd windA() const { return wind_->A(); }
	virtual vecd windB() const { return wind_->B(); }
	virtual vecd windC() const { return wind_->C(); }
// disk_irr_source_
protected:
	static std::shared_ptr<BasicFreddiIrradiationSource> initializeFreddiIrradiationSource(const std::string& angular_dist_type);
public:
	inline double angular_dist_disk(const double mu) const { return disk_irr_source_->angular_dist(mu); }
// opt_str_
protected:
	virtual void invalidate_optional_structure();

	template <DiskIntegrationRegion Region> size_t region_first() const {
		if constexpr(Region == HotRegion) {
			return first();
		} else if constexpr(Region == ColdRegion) {
			return last() + 1;
		} else {
			static_assert("Wrong Region template argument");
		}
	}
	template <DiskIntegrationRegion Region> size_t region_last() const {
		if constexpr(Region == HotRegion) {
			return last();
		} else if constexpr(Region == ColdRegion) {
			return Nx() - 1;
		} else {
			static_assert("Wrong Region template argument");
		}
	}
	template <DiskIntegrationRegion Region> double integrate(const vecd& values) const {
		return disk_radial_trapz(R(), values, region_first<Region>(), region_last<Region>());
	}
	template <DiskIntegrationRegion Region> double integrate(const std::function<double (size_t)>& func) const {
		return disk_radial_trapz(R(), func, region_first<Region>(), region_last<Region>());
	}
	template <DiskIntegrationRegion Region> double integrate(const vecd& x, const std::function<double (size_t)>& func) const {
		return trapz(x, func, region_first<Region>(), region_last<Region>());
	}
	template <DiskIntegrationRegion Region> double lazy_integrate(boost::optional<double>& x, const vecd& values) {
		if (!x) {
			x = integrate<Region>(values);
		}
		return *x;
	}
	template<DiskIntegrationRegion Region> double lazy_integrate(boost::optional<double> &opt, const vecd& x, const std::function<double (size_t)>& values) {
		if (!opt) {
			opt = integrate<Region>(x, values);
		}
		return *opt;
	}
	template <DiskIntegrationRegion Region> double I_lambda(double lambda) {
		const vecd* T;
		if constexpr(Region == HotRegion) {
			T = &Tph();
		} else if constexpr(Region == ColdRegion) {
			T = &Tirr();
		} else {
			static_assert("Wrong Region template argument");
		}
		return integrate<Region>([T, lambda](const size_t i) -> double { return Spectrum::Planck_lambda((*T)[i], lambda); });
	}
	double lazy_magnitude(boost::optional<double>& m, double lambda, double F0);
	virtual const vecd& Qx();
public:
	double Lx();
	const vecd& W();
	const vecd& Sigma();
	const vecd& Tph();
	const vecd& Tph_vis();
	const vecd& Tph_X();
	const vecd& Tirr();
	const vecd& Kirr();
	const vecd& Height();
	double Luminosity(const vecd& T, double nu1, double nu2) const;
	inline double magnitude(const double lambda, const double F0) {
		return -2.5 * std::log10(I_lambda<HotRegion>(lambda) * cosiOverD2() / F0);
	}
	template <DiskIntegrationRegion Region> double flux_region(double lambda) {
		return I_lambda<Region>(lambda) * m::pow<2>(lambda) / GSL_CONST_CGSM_SPEED_OF_LIGHT * cosiOverD2();
	}
	template <DiskIntegrationRegion Region> double flux_region(const Passband& passband) {
		const double intens = trapz(
				passband.lambdas,
				[this, &passband](const size_t i) -> double {
					return I_lambda<Region>(passband.lambdas[i]) * passband.transmissions[i];
				},
				0,
				passband.data.size() - 1);
		return intens * cosiOverD2() / passband.t_dnu;
	}
	inline double flux(const double lambda) { return flux_region<HotRegion>(lambda); }
	inline double flux(const Passband& passband) { return flux_region<HotRegion>(passband); }
	double flux_star(double lambda, double phase);
	double flux_star(const Passband& passband, double phase);
	inline double flux_star(double lambda) { return flux_star(lambda, phase_opt()); }
	inline double flux_star(const Passband& passband) { return flux_star(passband, phase_opt()); }
	inline double mU() { return lazy_magnitude(opt_str_.mU, lambdaU, irr0U); }
	inline double mB() { return lazy_magnitude(opt_str_.mB, lambdaB, irr0B); }
	inline double mV() { return lazy_magnitude(opt_str_.mV, lambdaV, irr0V); }
	inline double mR() { return lazy_magnitude(opt_str_.mR, lambdaR, irr0R); }
	inline double mI() { return lazy_magnitude(opt_str_.mI, lambdaI, irr0I); }
	inline double mJ() { return lazy_magnitude(opt_str_.mJ, lambdaJ, irr0J); }
	inline double Mdisk() { return lazy_integrate<HotRegion>(opt_str_.Mdisk, Sigma()); }
	double Mdot_wind();
	double Sigma_minus(double r) const;
	double Sigma_plus(double r) const;
	double R_cooling_front(double r);
	double v_cooling_front(double r);
};

#endif //FREDDI_FREDDI_STATE_HPP
