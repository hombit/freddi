# distutils: language = c++

from cython.operator cimport dereference

import numpy as np
cimport numpy as cnp

from freddi cimport *


cdef class Arguments:
    cdef FreddiArguments* cpp_args

    def __cinit__(
        self, *,
        double alpha=default_alpha, double Mx=default_Mx, double kerr=default_kerr, double Mopt=default_Mopt, double period=default_period, rin=None, rout=None,
        string opacity=default_opacity, string boundcond=default_boundcond, double Thot=default_Thot, string initialcond=default_initialcond, double F0=default_F0, double powerorder=default_powerorder, double gaussmu=default_gaussmu, double gausssigma=default_gausssigma, Mdisk0=None, Mdot0=None,
        double Cirr=default_Cirr, string irrfactortype=default_irrfactortype,
        double colourfactor=default_colourfactor, double emin=default_emin, double emax=default_emax, double inclination=default_inclination, double distance=default_distance, vector[double] lambdas=[],
        double time=default_time, double tau=default_tau, unsigned int Nx=default_Nx, string gridscale=default_gridscale, eps=None
    ):
        cdef GeneralArguments* general = new GeneralArguments(b'', b'', False)
        cdef BasicDiskBinaryArguments* basic
        if rin is None and rout is None:
            basic = new BasicDiskBinaryArguments(constructWithoutRinRout(alpha, Mx, kerr, Mopt, period))
        elif rin is None:
            basic = new BasicDiskBinaryArguments(constructWithoutRin(alpha, Mx, kerr, Mopt, period, rout))
        elif rout is None:
            basic = new BasicDiskBinaryArguments(constructWithoutRout(alpha, Mx, kerr, Mopt, period, rin))
        else:
            basic = new BasicDiskBinaryArguments(alpha, Mx, kerr, Mopt, period, rin, rout)
        cdef bint is_Mdisk0_specified = Mdisk0 is not None
        if Mdisk0 is None:
            Mdisk0 = -1.0
        cdef bint is_Mdot0_specified = Mdot0 is not None
        if Mdot0 is None:
            Mdot0 = -1.0
        cdef DiskStructureArguments* disk = new DiskStructureArguments(dereference(basic), opacity, boundcond, Thot, initialcond, F0, powerorder, gaussmu, gausssigma, is_Mdisk0_specified, is_Mdot0_specified, Mdisk0, Mdot0)
        cdef SelfIrradiationArguments* irr = new SelfIrradiationArguments(Cirr, irrfactortype)
        cdef FluxArguments* flux = new FluxArguments(colourfactor, emin, emax, inclination, distance, lambdas)
        cdef CalculationArguments* calc
        if eps is None:
            calc = new CalculationArguments(time, tau, Nx, gridscale)
        else:
            calc = new CalculationArguments(time, tau, Nx, gridscale, eps)
        self.cpp_args = new FreddiArguments(general, basic, disk, irr, flux, calc)

    def __dealloc__(self):
        del self.cpp_args

    @property
    def time(self):
        return self.cpp_args.calc.get().time


cdef class State:
    cdef FreddiState* cpp_state

    def __cinit__(self, Freddi freddi):
        self.cpp_state = new FreddiState(freddi.evolution.get_state())

    def __dealloc__(self):
        del self.cpp_state


cdef class Freddi:
    cdef FreddiEvolution* evolution
    cdef Arguments args

    def __cinit__(self, Arguments args):
        self.args = args
        self.evolution = new FreddiEvolution(dereference(self.args.cpp_args))

    def __dealloc__(self):
        del self.evolution

    cdef State get_state(self):
        return State(self)

    def __iter__(self):
        yield self.get_state()
        for tau in self.tau_view:
            self.evolution.step(tau)
            yield self.get_state()
            
