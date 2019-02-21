# distutils: language = c++
# cython: language_level = 3

from typing import Iterable
from cython.operator cimport dereference
from libc cimport math

import numpy as np
cimport numpy as cnp; cnp.import_array()

from freddi cimport *


cdef class State:
    """Disk radial structure

    Objects of this class shouldn't be created manually, normally they are
    obtained from `Freddi` methods

    """

    cdef FreddiState* cpp_state

    @property
    def Mdot_in(self) -> double:
       return self.cpp_state.Mdot_in()

    @property
    def Mdot(self) -> double:
        return self.Mdot_in

    @property
    def Mdot_out(self) -> double:
        return self.cpp_state.Mdot_out()

    @property
    def Lx(self) -> double:
        return self.cpp_state.Lx()
        
    @property
    def t(self) -> double:
        return self.cpp_state.t()

    @property
    def i_t(self) -> int:
        return self.cpp_state.i_t()

    @property
    def Nx(self) -> int:
        return self.cpp_state.Nx()

    @property
    def h(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.h().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def R(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.R().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def F(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.F().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def W(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.W().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tph(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Tph().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tph_vis(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Tph_vis().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tirr(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Tirr().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Cirr(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Cirr().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Sigma(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Sigma().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Height(self) -> np.ndarray[np.float]:
        cdef size_t first = self.cpp_state.first()
        cdef size_t last = self.cpp_state.last()
        cdef const double* data = self.cpp_state.Height().data() + first
        arr = np.asarray(<const double[:last + 1 - first]> data)
        arr.flags.writeable = False
        return arr

    @property
    def first(self) -> int:
        return self.cpp_state.first()

    @property
    def last(self) -> int:
        return self.cpp_state.last()

    @property
    def first_h(self) -> double:
        return self.cpp_state.h()[self.cpp_state.first()]

    @property
    def first_R(self) -> double:
        return self.cpp_state.R()[self.cpp_state.first()]

    @property
    def first_F(self) -> double:
        return self.cpp_state.F()[self.cpp_state.first()]

    @property
    def first_W(self) -> double:
        return self.cpp_state.W()[self.cpp_state.first()]

    @property
    def first_Tph(self) -> double:
        return self.cpp_state.Tph()[self.cpp_state.first()]

    @property
    def first_Tph_vis(self) -> double:
        return self.cpp_state.Tph_vis()[self.cpp_state.first()]

    @property
    def first_Tirr(self) -> double:
        return self.cpp_state.Tirr()[self.cpp_state.first()]

    @property
    def first_Cirr(self) -> double:
        return self.cpp_state.Cirr()[self.cpp_state.first()]

    @property
    def first_Sigma(self) -> double:
        return self.cpp_state.Sigma()[self.cpp_state.first()]

    @property
    def first_Height(self) -> double:
        return self.cpp_state.Height()[self.cpp_state.first()]

    @property
    def last_h(self) -> double:
        return self.cpp_state.h()[self.cpp_state.last()]

    @property
    def last_R(self) -> double:
        return self.cpp_state.R()[self.cpp_state.last()]

    @property
    def last_F(self) -> double:
        return self.cpp_state.F()[self.cpp_state.last()]

    @property
    def last_W(self) -> double:
        return self.cpp_state.W()[self.cpp_state.last()]

    @property
    def last_Tph(self) -> double:
        return self.cpp_state.Tph()[self.cpp_state.last()]

    @property
    def last_Tph_vis(self) -> double:
        return self.cpp_state.Tph_vis()[self.cpp_state.last()]

    @property
    def last_Tirr(self) -> double:
        return self.cpp_state.Tirr()[self.cpp_state.last()]

    @property
    def last_Cirr(self) -> double:
        return self.cpp_state.Cirr()[self.cpp_state.last()]

    @property
    def last_Sigma(self) -> double:
        return self.cpp_state.Sigma()[self.cpp_state.last()]

    @property
    def last_Height(self) -> double:
        return self.cpp_state.Height()[self.cpp_state.last()]

    @property
    def mU(self) -> double:
        return self.cpp_state.mU()

    @property
    def mB(self) -> double:
        return self.cpp_state.mB()

    @property
    def mV(self) -> double:
        return self.cpp_state.mV()

    @property
    def mR(self) -> double:
        return self.cpp_state.mR()

    @property
    def mI(self) -> double:
        return self.cpp_state.mI()

    @property
    def mJ(self) -> double:
        return self.cpp_state.mJ()

    @property
    def Mdisk(self) -> double:
        return self.cpp_state.Mdisk()

    def flux(self, lmbd) -> np.ndarray:
        arr = np.empty(np.broadcast(lmbd).shape, np.float)
        cdef cnp.broadcast it = cnp.broadcast(lmbd, arr)
        while cnp.PyArray_MultiIter_NOTDONE(it):
            (<double*>cnp.PyArray_MultiIter_DATA(it, 1))[0] = self.cpp_state.flux(
                (<double*>cnp.PyArray_MultiIter_DATA(it, 0))[0]
            )
            cnp.PyArray_MultiIter_NEXT(it)
        return arr
   

cdef class EvolutionResults:
    """Temporal distribution of various disk parameters

    Parameters
    ----------
    states: Iterable of states, e.g. list of Freddi
        States to be stored as results

    Methods
    -------
    flux(lmbd) : array
        Optical flux of the disk

    """

    cdef object[:] states

    def __cinit__(self, states: Iterable[State]):
        self.states = np.array(list(states), dtype=object)

    def flux(self, lmbd) -> np.ndarray:
        """Optical flux of the disk

        Parameters
        ----------
        lmbd : array_like
            Observation wavelength

        Returns
        -------
        array

        """
        cdef tuple lmbd_shape = cnp.broadcast(lmbd).shape
        arr = np.empty((self.states.size,) + cnp.broadcast(lmbd).shape, dtype=np.float)
        cdef size_t i
        for i in range(self.states.size):
            arr[i] = self.states[i].flux(lmbd)
        return arr

    def __getattr__(self, attr) -> np.ndarray:
        cdef size_t len_states = self.states.size
        first_val = np.asarray(getattr(self.states[0], attr))

        cdef size_t i

        if first_val.ndim == 0:
            arr = np.empty(len_states, dtype=first_val.dtype)
            for i in range(len_states):
                arr[i] = getattr(self.states[i], attr)
            return arr
        elif first_val.ndim == 1:
            arr = np.full((len_states, self.states[0].Nx), np.nan, dtype=first_val.dtype)
            for i in range(len_states):
                value = np.asarray(getattr(self.states[i], attr))
                arr[i, self.states[i].first:self.states[i].last+1] = value
            return arr
        else:
            raise ValueError("{}.ndim > 1 ({})".format(attr, first_val.ndim))


cdef class _BasicFreddi:
    """Accretion disk evolution modeller

    This object is iterable, it yields `State` object on each iteration step

    All parameters are keyword-only, all physical values are in CGS units.
    See detailed description of parameters in documentation

    Parameters
    ----------
    alpha : float, optional
    Mx : float, optional
    kerr : float, optional
    Mopt : float, optional
    period : float, optional
    rin : float or None, optional
    rout : float or None, optional
    opacity : bytes, optional
        Should be `b'Kramers'` or `b'OPAL'`
    boundcond : bytes, optional
        Should be `b'Teff'` or `b'Tirr'`
    Mdotout : float, optional
    Thot : float, optional
    initialcond : bytes, optional
        Should be `b'powerF'`, `b'powerSigma'`, `b'sinusF'`, `b'gaussF'` or
        `b'quasistat'`
    FO : float, optional
    Mdisk0 : float, optional
    Mdot0 : float, optional
    powerorder : float, optional
    gaussmu : float, optional
    gausssigma : float, optional
    wind : bytes, optional
    windparams : list of float, optional
    Cirr : float, optional
    irrfactortype : bytes, optional
        Should be `b'const'` or `b'square'`
    colourfactor : float, optional
    emin : float, optional
    emax : float, optional
    inclination : float, optional
        In degrees
    distance : float, optional
    lambda : list of float, optional
    time : float, optional
    tau : float, optional
    Nx : int, optional
    gridscale : bytes, optional
        Should be `b'log'` or `b'linear'`{add_args}

    Attributes
    ----------
    time : float
        Specified time
    distance : float
        Specified distance
    alpha : float
        Shakura-Sunyaev turbulent parameter
    cosi : flaot
        Cosinus of specified inclination
    lambdas : numpy.ndarray
        Specified lambdas, in cm
    Nt : int
        Number of evolutions steps. The number of time moments where disk
        structure will be obtained is larger by unity
    Cirr : float
        Irradiation factor. Can be changed via assignment operator

    Methods
    -------
    evolve() : EvolutionResults
        Calculate disk evolution
    alt(**kwargs) : Freddi
        Alternative Freddi constructor accepting astropy Quantity and str

    """

    cdef FreddiArguments* args
    cdef FreddiEvolution* evolution

    def __cinit__(
        self, *, bint cgs=True, bytes type=b'Normal',
        double alpha=default_alpha, double Mx=default_Mx, double kerr=default_kerr, double Mopt=default_Mopt,
            double period=default_period, rin=None, rout=None,
        string opacity=default_opacity, double Mdotout=default_Mdotout,
            string boundcond=default_boundcond, double Thot=default_Thot, string initialcond=default_initialcond,
            double F0=default_F0, double powerorder=default_powerorder, double gaussmu=default_gaussmu,
            double gausssigma=default_gausssigma, Mdisk0=None, Mdot0=None,
            string wind=default_wind, vector[double] windparams=[],
        double Cirr=default_Cirr, string irrfactortype=default_irrfactortype,
        double colourfactor=default_colourfactor, double emin=default_emin, double emax=default_emax,
            double inclination=default_inclination, double distance=default_distance, vector[double] lambdas=[],
        double time=default_time, double tau=default_tau, unsigned int Nx=default_Nx,
            string gridscale=default_gridscale, eps=None,
        **kwargs,
    ):
        if not cgs:
            Mx = sunToGram(Mx)
            Mopt = sunToGram(Mopt)
            period = dayToS(period)
            if rin is not None:
                rin = rgToCm(rin, Mx)
            if rout is not None:
                rout = sunToCm(rout)
            emin = kevToHertz(emin)
            emax = kevToHertz(emax)
            distance = kpcToCm(distance)
            for i in range(lambdas.size()):
                lambdas[i] = angstromToCm(lambdas[i])
            time = dayToS(time)
            tau = dayToS(tau)
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
        cdef DiskStructureArguments* disk = new DiskStructureArguments(
            dereference(basic), opacity, Mdotout, boundcond, Thot, initialcond, F0, powerorder,
            gaussmu, gausssigma, is_Mdisk0_specified, is_Mdot0_specified, Mdisk0, Mdot0, wind, windparams
        )
        cdef SelfIrradiationArguments* irr = new SelfIrradiationArguments(Cirr, irrfactortype)
        cdef FluxArguments* flux = new FluxArguments(colourfactor, emin, emax, inclination, distance, lambdas)
        cdef CalculationArguments* calc
        if eps is None:
            calc = new CalculationArguments(time, tau, Nx, gridscale)
        else:
            calc = new CalculationArguments(time, tau, Nx, gridscale, eps)
        self.args = new FreddiArguments(general, basic, disk, irr, flux, calc)

    def __dealloc__(self):
        del self.args
        del self.evolution

    @classmethod
    def alt(cls, **kwargs):
        """Alternative Freddi constructor accepting astropy Quantity and str

        All physical values can be passed as `astropy.units.Quantity`, and all
        string values can be passed as `str`. Needs `astropy` module to run

        Parameters
        ----------
        **kwargs :
            `Freddi` arguments

        Returns
        -------
        Freddi

        """
        from astropy.units import Quantity

        for item, value in kwargs.items():
            if isinstance(value, Quantity):
                kwargs[item] = value.cgs.value
                continue
            if isinstance(value, str):
                kwargs[item] = value.encode()
        return cls(**kwargs)

    @property
    def time(self) -> double:
        return self.args.calc.get().time

    @property
    def distance(self) -> double:
        return self.args.flux.get().distance

    @property
    def alpha(self) -> double:
        return self.args.basic.get().alpha

    @property
    def cosi(self) -> double:
        return math.cos(self.args.flux.get().inclination / (<double> 180.) * math.M_PI)

    @property
    def lambdas(self) -> np.ndarray:
        cdef const double* data = self.args.flux.get().lambdas.data()
        cdef size_t size = self.args.flux.get().lambdas.size()
        if size == (<size_t> 0):
            return np.array([], dtype=np.float)
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Nt(self) -> int:
        return self.evolution.state().Nt()

    cdef void change_SelfIrradiationArguments(self, Cirr=None, irrfactortype=None):
        cdef double c_Cirr = self.args.irr.get().Cirr if Cirr is None else Cirr
        cdef string c_irrfactortype = self.args.irr.get().irrfactortype if irrfactortype is None else irrfactortype
        cdef SelfIrradiationArguments* irr = new SelfIrradiationArguments(c_Cirr, c_irrfactortype)
        self.args.irr.reset(irr)

    @property
    def Cirr(self) -> double:
        return self.args.irr.get().Cirr

    @Cirr.setter
    def Cirr(self, val: double) -> None:
        self.change_SelfIrradiationArguments(Cirr=val)

    cdef void change_DiskStructureArguments(self, opacity=None, Mdotout=None, boundcond=None, Thot=None):
            cdef string c_opacity = self.args.disk.get().opacity if opacity is None else opacity
            cdef double c_Mdotout = self.args.disk.get().Mdotout if Mdotout is None else Mdotout
            cdef string c_boundcond = self.args.disk.get().boundcond if boundcond is None else boundcond
            cdef double c_Thot = self.args.disk.get().Thot if Thot is None else Thot
            cdef DiskStructureArguments* disk = new DiskStructureArguments(
                dereference(self.args.basic.get()),
                c_opacity, c_Mdotout, c_boundcond, c_Thot,
                self.args.disk.get().initialcond, self.args.disk.get().F0, self.args.disk.get().powerorder,
                self.args.disk.get().gaussmu, self.args.disk.get().gausssigma,
                True, True, self.args.disk.get().Mdisk0, self.args.disk.get().Mdot0,
                self.args.disk.get().wind, self.args.disk.get().windparams,
            )
            self.args.disk.reset(disk)

    @property
    def boundcond(self) -> bytes:
        return self.args.disk.get().boundcond

    @boundcond.setter
    def boundcond(self, val: bytes) -> None:
        self.change_DiskStructureArguments(opacity=None, Mdotout=None, boundcond=val)

    @property
    def Thot(self) -> double:
        return self.args.disk.get().Thot

    @Thot.setter
    def Thot(self, val: double) -> None:
        self.change_DiskStructureArguments(opacity=None, Mdotout=None, boundcond=None, Thot=val)

    def get_state(self):
        raise NotImplementedError

    def __iter__(self) -> State:
        """Iterate disk over time

        The first yielded object is for initial distribution

        Yields
        ------
        State
            `State` disk radial structure object

        """
        cdef State state = self.get_state()
        while state.t <= self.time:
            yield state
            self.evolution.step()
            state = self.get_state()

    def evolve(self):
        """Calculate disk evolution

        Returns
        -------
        EvolutionResults

        """
        return EvolutionResults(self)

    def __getattr__(self, attr):
        state = self.get_state()
        return getattr(state, attr)


cdef class Freddi(_BasicFreddi):
    __doc__ = _BasicFreddi.__doc__.format(add_args='')

    def __cinit__(self, **kwargs):
        self.evolution = new FreddiEvolution(dereference(self.args))

    cpdef State get_state(self):
        cdef State state = State.__new__(State)
        state.cpp_state = <FreddiState*> new FreddiEvolution(dereference(self.evolution))
        return state


cdef class FreddiNeutronStar(_BasicFreddi):
    __doc__ = _BasicFreddi.__doc__.format(add_args="""
    Rx : float, optional
    freqx : float, optional
    Bx : float, optional
    epsilonAlfven : float, optional""")

    cdef FreddiNeutronStarEvolution* ns_evolution

    def __cinit__(
        self, *,
        double Rx=default_Rx, double freqx=default_freqx, double Bx=default_Bx, double epsilonAlfven=default_epsilonAlfven, double inversebeta=default_inversebeta, double Rdead=default_Rdead,
        **kwargs,
    ):
        cdef NeutronStarArguments* ns = new NeutronStarArguments(Rx, freqx, Bx, epsilonAlfven, inversebeta, Rdead)
        cdef FreddiNeutronStarArguments* ns_args = new FreddiNeutronStarArguments(dereference(self.args), ns)
        self.args = <FreddiArguments*> ns_args
        self.ns_evolution = new FreddiNeutronStarEvolution(dereference(ns_args))
        self.evolution = <FreddiEvolution*> self.ns_evolution

    cpdef State get_state(self, bint shadow=False):
        cdef State state = State.__new__(State)
        state.cpp_state = <FreddiState*> new FreddiNeutronStarEvolution(dereference(self.ns_evolution))
        return state

    @property
    def Fmagn(self) -> np.ndarray[np.float]:
        cdef const double* data = self.ns_evolution.Fmagn.data()
        cdef size_t size = self.ns_evolution.Fmagn.size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def dFmagn_dh(self) -> np.ndarray[np.float]:
        cdef const double* data = self.ns_evolution.dFmagn_dh.data()
        cdef size_t size = self.ns_evolution.dFmagn_dh.size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def d2Fmagn_dh2(self) -> np.ndarray[np.float]:
        cdef const double* data = self.ns_evolution.d2Fmagn_dh2.data()
        cdef size_t size = self.ns_evolution.d2Fmagn_dh2.size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def R_cor(self) -> double:
        return self.ns_evolution.R_cor
