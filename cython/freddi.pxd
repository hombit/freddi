# cython: language_level = 3

from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from 'freddi_state.hpp':
    cdef cppclass FreddiState:
        FreddiState(const FreddiState&)
        size_t Nt
        size_t Nx
        size_t first()
        size_t last()
        double Mdot_in()
        double Mdot_out()
        double Lx()
        double t()
        size_t i_t()
        const vector[double]& h()
        const vector[double]& R()
        const vector[double]& F()
        const vector[double]& W()
        const vector[double]& Tph()
        const vector[double]& Tph_vis()
        const vector[double]& Tirr()
        const vector[double]& Cirr()
        const vector[double]& Sigma()
        const vector[double]& Height()
        double flux(double)
        double mU() 
        double mB() 
        double mV() 
        double mR() 
        double mI() 
        double mJ() 
        double Mdisk()

cdef extern from 'freddi_evolution.hpp':
    cdef cppclass FreddiEvolution:
        FreddiEvolution(const FreddiArguments&) except +
        FreddiEvolution(const FreddiEvolution&)
        void step() except +
        void step(double) except +
        FreddiState& state()
    cdef cppclass FreddiNeutronStarEvolution:
        FreddiNeutronStarEvolution(const FreddiNeutronStarArguments&) except +
        FreddiNeutronStarEvolution(const FreddiNeutronStarEvolution&)
        const vector[double] Fmagn
        const vector[double] dFmagn_dh
        const vector[double] d2Fmagn_dh2
        const double R_cor


cdef extern from 'arguments.hpp':
    cdef cppclass GeneralArguments:
        GeneralArguments(const string&, const string&, bint)
    cdef cppclass BasicDiskBinaryArguments:
        BasicDiskBinaryArguments(double, double, double, double, double, double, double) except +
        BasicDiskBinaryArguments(BasicDiskBinaryArguments&&)
        const double alpha, Mx, kerr, Mopt, period, rin, rout
    cdef cppclass DiskStructureArguments:
        DiskStructureArguments(const BasicDiskBinaryArguments&, const string&, double, const string&, double, const string&, double, double, double, double, bint, bint, double, double, const string&, const vector[double]&) except +
        const string opacity, boundcond, initialcond, wind
        const double Mdotout, Thot, F0, powerorder, gaussmu, gausssigma, Mdisk0, Mdot0
        const vector[double] windparams
    cdef cppclass SelfIrradiationArguments:
        SelfIrradiationArguments(double, const string&) except +
        const double Cirr
        const string irrfactortype
    cdef cppclass FluxArguments:
        FluxArguments(double, double, double, double, double, vector[double])
        const double colourfactor, emin, emax, inclination, distance
        const vector[double] lambdas
    cdef cppclass CalculationArguments:
        CalculationArguments(double, double, unsigned int, const string&) except +
        CalculationArguments(double, double, unsigned int, const string&, double) except +
        const double time, tau, eps
        const unsigned int Nx
        const string gridscale
    cdef cppclass FreddiArguments:
        FreddiArguments(GeneralArguments* general, BasicDiskBinaryArguments* basic, DiskStructureArguments* disk, SelfIrradiationArguments* irr, FluxArguments* flux, CalculationArguments* cal)
        shared_ptr[GeneralArguments] general
        shared_ptr[BasicDiskBinaryArguments] basic
        shared_ptr[DiskStructureArguments] disk
        shared_ptr[SelfIrradiationArguments] irr
        shared_ptr[FluxArguments] flux
        shared_ptr[CalculationArguments] calc
    cdef cppclass NeutronStarArguments:
        NeutronStarArguments(double, double, double, double, double, double)
        const double Rx, freqx, Bx, epsilonAlfven, Fdead
    cdef cppclass FreddiNeutronStarArguments:
        FreddiNeutronStarArguments(const FreddiArguments& freddi_args, NeutronStarArguments* ns)
        shared_ptr[NeutronStarArguments] ns


cdef extern from 'arguments.hpp' namespace 'BasicDiskBinaryArguments':
    cdef BasicDiskBinaryArguments constructWithoutRinRout(double, double, double, double, double)
    cdef BasicDiskBinaryArguments constructWithoutRin(double, double, double, double, double, double)
    cdef BasicDiskBinaryArguments constructWithoutRout(double, double, double, double, double, double)
    cdef const double default_alpha
    cdef const double default_Mx
    cdef const double default_kerr
    cdef const double default_accfreq
    cdef const double default_Mopt
    cdef const double default_period


cdef extern from 'arguments.hpp' namespace 'DiskStructureArguments':
    cdef const char* default_opacity
    cdef const double default_Mdotout
    cdef const char* default_boundcond
    cdef const double default_Thot
    cdef const char* default_initialcond
    cdef const double default_F0
    cdef const double default_powerorder
    cdef const double default_gaussmu
    cdef const double default_gausssigma
    cdef const char* default_wind


cdef extern from 'arguments.hpp' namespace 'SelfIrradiationArguments':
    cdef const double default_Cirr
    cdef const char* default_irrfactortype


cdef extern from 'arguments.hpp' namespace 'FluxArguments':
    cdef const double default_colourfactor
    cdef const double default_emin
    cdef const double default_emax
    cdef const double default_inclination
    cdef const double default_distance


cdef extern from 'arguments.hpp' namespace 'CalculationArguments':
    cdef const double default_time
    cdef const double default_tau
    cdef const unsigned int default_Nx
    cdef const char* default_gridscale


cdef extern from 'arguments.hpp' namespace 'NeutronStarArguments':
    cdef const double default_Rx, default_freqx, default_Bx, default_epsilonAlfven, default_inversebeta, default_Fdead


cdef extern from 'unit_transformation.hpp':
    cdef double sunToGram(double)
    cdef double rgToCm(double, double)
    cdef double sunToCm(double)
    cdef double kpcToCm(double)
    cdef double angstromToCm(double)
    cdef double dayToS(double)
    cdef double kevToHertz(double)
