from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from 'freddi.hpp':
    cdef cppclass FreddiState:
        FreddiState(const FreddiState&)
        double get_Mdot_in()
        double get_Mdot_out()
        double get_Lx()
        double get_t()
        size_t get_i_t()
        unsigned int get_Nx()
        const vector[double]& get_h() 
        const vector[double]& get_R() 
        const vector[double]& get_F() 
        const vector[double]& get_W() 
        const vector[double]& get_Tph() 
        const vector[double]& get_Tph_vis() 
        const vector[double]& get_Tirr()
        const vector[double]& get_Cirr() 
        const vector[double]& get_Sigma() 
        const vector[double]& get_Height() 
        double flux(double)
        double mU() 
        double mB() 
        double mV() 
        double mR() 
        double mI() 
        double mJ() 
        double Mdisk()
    cdef cppclass FreddiEvolution:
        FreddiEvolution(const FreddiArguments&) except +
        void step() except +
        void step(double) except +
        vector[FreddiState] evolve() except +
        const FreddiState& get_state()
        size_t Nt


cdef extern from 'arguments.hpp':
    cdef cppclass GeneralArguments:
        GeneralArguments(const string&, const string&, bint)
    cdef cppclass BasicDiskBinaryArguments:
        BasicDiskBinaryArguments(double, double, double, double, double, double, double) except +
        BasicDiskBinaryArguments(BasicDiskBinaryArguments&&)
        const double alpha, Mx, kerr, Mopt, period, rin, rout
    cdef cppclass DiskStructureArguments:
        DiskStructureArguments(const BasicDiskBinaryArguments&, const string&, const string&, double, const string&, double, double, double, double, bint, bint, double, double) except +
        const string opacity, boundcond, initialcond
        const double Thot, F0, powerorder, gaussmu, gausssigma, Mdisk0, Mdot0
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


cdef extern from 'arguments.hpp' namespace 'BasicDiskBinaryArguments':
    cdef BasicDiskBinaryArguments constructWithoutRinRout(double, double, double, double, double)
    cdef BasicDiskBinaryArguments constructWithoutRin(double, double, double, double, double, double)
    cdef BasicDiskBinaryArguments constructWithoutRout(double, double, double, double, double, double)
    cdef const double default_alpha
    cdef const double default_Mx
    cdef const double default_kerr
    cdef const double default_Mopt
    cdef const double default_period


cdef extern from 'arguments.hpp' namespace 'DiskStructureArguments':
    cdef const char* default_opacity
    cdef const char* default_boundcond
    cdef const double default_Thot
    cdef const char* default_initialcond
    cdef const double default_F0
    cdef const double default_powerorder
    cdef const double default_gaussmu
    cdef const double default_gausssigma


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


cdef extern from 'unit_transformation.hpp':
    cdef double sunToGram(double)
    cdef double rgToCm(double, double)
    cdef double sunToCm(double)
    cdef double kpcToCm(double)
    cdef double angstromToCm(double)
    cdef double dayToS(double)
    cdef double kevToHertz(double)
