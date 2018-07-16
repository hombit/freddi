# distutils: language = c++

from cython.operator cimport dereference

from freddi cimport *

cdef class Arguments:
    cdef GeneralArguments* general
    cdef BasicDiskBinaryArguments* basic
    cdef DiskStructureArguments* disk
    cdef SelfIrradiationArguments* irr
    cdef FluxArguments* flux
    cdef CalculationArguments* calc
    cdef FreddiArguments* cpp_args

    def __cinit__(
        self, *,
        double alpha=default_alpha, double Mx=default_Mx, double kerr=default_kerr, double Mopt=default_Mopt, double period=default_period, rin=None, rout=None,
        string opacity=default_opacity, string boundcond=default_boundcond, double Thot=default_Thot, string initialcond=default_initialcond, double F0=default_F0, double powerorder=default_powerorder, double gaussmu=default_gaussmu, double gausssigma=default_gausssigma, Mdisk0=None, Mdot0=None,
        double Cirr=default_Cirr, string irrfactortype=default_irrfactortype,
        double colourfactor=default_colourfactor, double emin=default_emin, double emax=default_emax, double inclination=default_inclination, double distance=default_distance, vector[double] lambdas=[],
        double time=default_time, double tau=default_tau, unsigned int Nx=default_Nx, string gridscale=default_gridscale, eps=None
    ):
        self.general = new GeneralArguments(b'', b'', False)
        if rin is None and rout is None:
            self.basic = new BasicDiskBinaryArguments(constructWithoutRinRout(alpha, Mx, kerr, Mopt, period))
        elif rin is None:
            self.basic = new BasicDiskBinaryArguments(constructWithoutRin(alpha, Mx, kerr, Mopt, period, rout))
        elif rout is None:
            self.basic = new BasicDiskBinaryArguments(constructWithoutRout(alpha, Mx, kerr, Mopt, period, rin))
        else:
            self.basic = new BasicDiskBinaryArguments(alpha, Mx, kerr, Mopt, period, rin, rout)
        cdef bint is_Mdisk0_specified = Mdisk0 is not None
        if Mdisk0 is None:
            Mdisk0 = -1.0
        cdef bint is_Mdot0_specified = Mdot0 is not None
        if Mdot0 is None:
            Mdot0 = -1.0
        self.disk = new DiskStructureArguments(dereference(self.basic), opacity, boundcond, Thot, initialcond, F0, powerorder, gaussmu, gausssigma, is_Mdisk0_specified, is_Mdot0_specified, Mdisk0, Mdot0)
        self.irr = new SelfIrradiationArguments(Cirr, irrfactortype)
        self.flux = new FluxArguments(colourfactor, emin, emax, inclination, distance, lambdas)
        if eps is None:
            self.calc = new CalculationArguments(time, tau, Nx, gridscale)
        else:
            self.calc = new CalculationArguments(time, tau, Nx, gridscale, eps)
        self.cpp_args = new FreddiArguments(self.general, self.basic, self.disk, self.irr, self.flux, self.calc)

    @property
    def time(self):
        return self.cpp_args.calc.get().time

