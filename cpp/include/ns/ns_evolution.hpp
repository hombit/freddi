#ifndef FREDDI_NS_EVOLUTION_HPP
#define FREDDI_NS_EVOLUTION_HPP

#include <freddi_evolution.hpp>
#include <ns/ns_arguments.hpp>

class FreddiNeutronStarEvolution: public FreddiEvolution {
private:
	class NeutronStarStructure {
	public:
		NeutronStarArguments args_ns;
		double k_t = 1. / 3.;
//      Factor previously used changing initial radius of the disk
//		double xi = 1.0;
//		double xi_pow_minus_7_2;
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
		virtual double operator()(double R_to_Rcor) const = 0;
	};

	class NoOutflowNSMdotFraction: public BasicNSMdotFraction {
	public:
		NoOutflowNSMdotFraction() = default;
		~NoOutflowNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	class PropellerNSMdotFraction: public BasicNSMdotFraction {
	public:
		PropellerNSMdotFraction() = default;
		~PropellerNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	class CorotationBlockNSMdotFraction: public BasicNSMdotFraction {
	public:
		CorotationBlockNSMdotFraction() = default;
		~CorotationBlockNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	// https://arxiv.org/pdf/1010.1528.pdf Eksi-Kutlu (2010)
	class EksiKultu2010NSMdotFraction: public BasicNSMdotFraction {
	public:
		EksiKultu2010NSMdotFraction() = default;
		~EksiKultu2010NSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	class Romanova2018NSMdotFraction: public BasicNSMdotFraction {
	private:
		const double par1;
		const double par2;
	public:
		Romanova2018NSMdotFraction(double par1, double par2);
		~Romanova2018NSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	class GeometricalNSMdotFraction: public BasicNSMdotFraction {
	private:
		const double chi;
	public:
		GeometricalNSMdotFraction(double chi);
		~GeometricalNSMdotFraction() override = default;
		virtual double operator()(double R_to_Rcor) const override;
	};

	class BasicNSAccretionEfficiency {
	public:
		virtual ~BasicNSAccretionEfficiency() = 0;
		virtual double operator()(const FreddiNeutronStarEvolution& freddi, double Rm) const;
		virtual double RmIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const = 0;
		virtual double RxIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const = 0;
		virtual double RiscoIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const = 0;
	};

	class DummyNSAccretionEfficiency: public BasicNSAccretionEfficiency {
	protected:
		double newtonian(const FreddiNeutronStarEvolution& freddi, double Rm) const;
	public:
		~DummyNSAccretionEfficiency() override = default;
		double RmIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return newtonian(freddi, Rm); }
		double RxIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return newtonian(freddi, Rm); }
		double RiscoIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return newtonian(freddi, Rm); }
	};

	class SibgatullinSunyaev2000NSAccretionEfficiency: public BasicNSAccretionEfficiency {
	protected:
		double schwarzschild(const FreddiNeutronStarEvolution& freddi, double Rm) const;
		double small_magnetosphere(const FreddiNeutronStarEvolution& freddi, double Rm) const;
	public:
		using BasicNSAccretionEfficiency::BasicNSAccretionEfficiency;
		~SibgatullinSunyaev2000NSAccretionEfficiency() override = default;
		double RmIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return schwarzschild(freddi, Rm); }
		double RxIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return small_magnetosphere(freddi, Rm); }
		double RiscoIsFurthest(const FreddiNeutronStarEvolution& freddi, double Rm) const override { return small_magnetosphere(freddi, Rm); }
	};
private:
	std::shared_ptr<const NeutronStarStructure> ns_str_;
	NeutronStarOptionalStructure ns_opt_str_;
	std::shared_ptr<BasicNSMdotFraction> fp_;
	std::shared_ptr<BasicNSAccretionEfficiency> eta_ns_;
private:
	static std::shared_ptr<BasicNSMdotFraction> initializeNsMdotFraction(const NeutronStarArguments& args_ns);
	static std::shared_ptr<BasicNSAccretionEfficiency> initializeNsAccretionEfficiency(const NeutronStarArguments& args_ns, const FreddiNeutronStarEvolution* freddi);
// ns_str_
public:
	inline double k_t() const { return ns_str_->k_t; }
//	inline double xi() const { return ns_str_->xi; }
//	inline double xi_pow_minus_7_2() const { return ns_str_->xi_pow_minus_7_2; }
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
	inline double fp(double radius) const { return (*fp_)(radius / R_cor()); }
	inline double fp() const { return fp(R()[first()]); }
// eta_ns_
public:
	inline double eta_ns(double Rm) const { return (*eta_ns_)(*this, Rm); }
	inline double eta_ns() const { return eta_ns(R_alfven()); }
public:
	using FreddiEvolution::step;
	virtual double Mdot_in() const override;
	double R_alfven() const;
protected:
	virtual void invalidate_optional_structure() override;
	virtual void truncateInnerRadius() override;
	virtual vecd windC() const override;
	virtual IrradiatedStar::sources_t star_irr_sources() override;
	virtual const vecd& Tph_X() override;
public:
	FreddiNeutronStarEvolution(const FreddiNeutronStarArguments& args);
	explicit FreddiNeutronStarEvolution(const FreddiNeutronStarEvolution&) = default;
	virtual const vecd& Qx() override;
	virtual double Lbol_disk() const override;
public:
	using iterator = EvolutionIterator<FreddiNeutronStarEvolution>;
	inline iterator begin() { return {this}; }
	inline iterator end() { return {Nt() + 1}; }
};

#endif //FREDDI_NS_EVOLUTION_HPP
