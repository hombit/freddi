# distutils: language = c++

from cython.operator cimport dereference

import numpy as np
cimport numpy as cnp; cnp.import_array()

from freddi cimport *


cdef class State:
    """Disk radial structure

    Objects of this class shouldn't be created manually, normally they are
    obtained from `Freddi` methods

    """

    cdef FreddiState* cpp_state

    def __dealloc__(self):
        if self.cpp_state:
            del self.cpp_state

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
        cdef const double* data = self.cpp_state.h().data()
        cdef size_t size = self.cpp_state.h().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def R(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.R().data()
        cdef size_t size = self.cpp_state.R().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def F(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.F().data()
        cdef size_t size = self.cpp_state.F().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def W(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.W().data()
        cdef size_t size = self.cpp_state.W().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tph(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Tph().data()
        cdef size_t size = self.cpp_state.Tph().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tph_vis(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Tph_vis().data()
        cdef size_t size = self.cpp_state.Tph_vis().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Tirr(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Tirr().data()
        cdef size_t size = self.cpp_state.Tirr().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Cirr(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Cirr().data()
        cdef size_t size = self.cpp_state.Cirr().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Sigma(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Sigma().data()
        cdef size_t size = self.cpp_state.Sigma().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def Height(self) -> np.ndarray[np.float]:
        cdef const double* data = self.cpp_state.Height().data()
        cdef size_t size = self.cpp_state.Height().size()
        arr = np.asarray(<const double[:size]> data)
        arr.flags.writeable = False
        return arr

    @property
    def last_h(self) -> double:
        return self.cpp_state.h().back()

    @property
    def last_R(self) -> double:
        return self.cpp_state.R().back()

    @property
    def last_F(self) -> double:
        return self.cpp_state.F().back()

    @property
    def last_W(self) -> double:
        return self.cpp_state.W().back()

    @property
    def last_Tph(self) -> double:
        return self.cpp_state.Tph().back()

    @property
    def last_Tph_vis(self) -> double:
        return self.cpp_state.Tph_vis().back()

    @property
    def last_Tirr(self) -> double:
        return self.cpp_state.Tirr().back()

    @property
    def last_Cirr(self) -> double:
        return self.cpp_state.Cirr().back()

    @property
    def last_Sigma(self) -> double:
        return self.cpp_state.Sigma().back()

    @property
    def last_Height(self) -> double:
        return self.cpp_state.Height().back()

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
   

cdef State state_from_cpp(const FreddiState& cpp_state):
    """Construct State object from C++ FreddiState"""
    cdef State state = State()
    state.cpp_state = new FreddiState(cpp_state)
    return state


cdef class EvolutionResults:
    """Temporal distribution of various disk parameters

    Parameters
    ----------
    freddi: Freddi
        Freddi object to be evolved and prepared

    Methods
    -------
    flux(lmbd) : array
        Optical flux of the disk

    """

    cdef vector[FreddiState] cpp_states
    cdef object[:] states

    def __cinit__(self, Freddi freddi):
        self.cpp_states = freddi.evolution.evolve()
        self.states = np.empty(<Py_ssize_t> self.cpp_states.size(), dtype=object)
        cdef Py_ssize_t i
        for i in range(self.states.size):
            self.states[i] = state_from_cpp(self.cpp_states[i])

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
        cdef size_t len_states = self.cpp_states.size()
        first_val = np.asarray(getattr(self.states[0], attr))
        cdef tuple shape_val = first_val.shape
        cdef tuple shape = (len_states,) + shape_val
        arr = np.full(shape, np.nan, dtype=first_val.dtype)
        arr_flat = arr.reshape(-1)
        cdef size_t i
        cdef size_t size_arr = arr.size
        cdef size_t len_line = size_arr / len_states
        for i in range(len_states):
            value = np.asarray(getattr(self.states[i], attr))
            arr_flat[i * len_line : i * len_line + value.size] = value.reshape(-1)
        return arr


cdef class Freddi:
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
    Mdotout: float, optional
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
        Should be `b'log'` or `b'linear'`

    Attributes
    ----------
    time : float
        Specified time
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
        self, *, bint cgs=True,
        double alpha=default_alpha, double Mx=default_Mx, double kerr=default_kerr, double Mopt=default_Mopt, double period=default_period, rin=None, rout=None,
        string opacity=default_opacity, double Mdotout=default_Mdotout, string boundcond=default_boundcond, double Thot=default_Thot, string initialcond=default_initialcond, double F0=default_F0, double powerorder=default_powerorder, double gaussmu=default_gaussmu, double gausssigma=default_gausssigma, Mdisk0=None, Mdot0=None,
        double Cirr=default_Cirr, string irrfactortype=default_irrfactortype,
        double colourfactor=default_colourfactor, double emin=default_emin, double emax=default_emax, double inclination=default_inclination, double distance=default_distance, vector[double] lambdas=[],
        double time=default_time, double tau=default_tau, unsigned int Nx=default_Nx, string gridscale=default_gridscale, eps=None
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
        cdef DiskStructureArguments* disk = new DiskStructureArguments(dereference(basic), opacity, Mdotout, boundcond, Thot, initialcond, F0, powerorder, gaussmu, gausssigma, is_Mdisk0_specified, is_Mdot0_specified, Mdisk0, Mdot0)
        cdef SelfIrradiationArguments* irr = new SelfIrradiationArguments(Cirr, irrfactortype)
        cdef FluxArguments* flux = new FluxArguments(colourfactor, emin, emax, inclination, distance, lambdas)
        cdef CalculationArguments* calc
        if eps is None:
            calc = new CalculationArguments(time, tau, Nx, gridscale)
        else:
            calc = new CalculationArguments(time, tau, Nx, gridscale, eps)
        self.args = new FreddiArguments(general, basic, disk, irr, flux, calc)

        self.evolution = new FreddiEvolution(dereference(self.args))

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
        return self.evolution.Nt

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

    cdef State get_state(self):
        return state_from_cpp(self.evolution.state())

    def __iter__(self):
        """Iterate disk over time

        The first yielded object is for initial distribution

        Yields
        ------
        State
            `State` disk radial structure object

        """
        state = self.get_state()
        yield state
        while state.t <= self.time:
            self.evolution.step()
            state = self.get_state()
            yield state

    def evolve(self):
        """Calculate disk evolution

        Returns
        -------
        EvolutionResults

        """
        return EvolutionResults(self)
