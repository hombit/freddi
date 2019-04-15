import unittest

import numpy as np
from numpy.testing import assert_allclose

from freddi import Freddi

try:
    import scipy.special
except ImportError:
    scipy = None


DAY = 86400


class ShakuraSunyaevSubctriticalTestCase(unittest.TestCase):
    def test(self):
        """Shakura-Sunyaev 1973 stationary subcritical disk

        The analytical stationary solution is $F ~ (h - h_{in})$

        The initial condition is non-stationary and does not satisfies outer
        boundary condition.

        """
        Mdot = 1e18
        fr = Freddi(initialcond=b'sinusF', Mdot0=Mdot/2, Mdotout=Mdot, time=1000*DAY, tau=1*DAY)
        for state in fr:
            pass
        h = state.h
        assert_allclose(state.F, Mdot * (h - h[0]), rtol=1e-5)


class ShakuraSunyaevSupercriticalTestCase(unittest.TestCase):
    def test(self):
        """Shakura-Sunayev 1973 stationary supercritical disk with outflow

        $d\dot M_w / dA = -\dot M_{crit} / (2 \pi R_{sph} r)$

        The analytical stationary solution is $F ~ (h^3 - h^3_{in})$ for
        the region where disk is supercritical: radius $r$ is between $R_{in}$
        and $R_{sph} = R_{in} \dot M_{out} / \dot M_{crit}$, $Mdot_{crit}$ is
        critical (Eddington) accretion rate on the black hole.

        The initial condition is non-stationary and does not satisfies outer
        boundary condition, outer radius is $R_{sph}$ for given accretion rate
        through this radius.

        """
        Mx = 2e34
        GM = 6.673e-8 * Mx
        c = 2.99792458e10
        Rin = 6 * 6.673e-8 * Mx / c**2
        Rout = 100 * Rin
        Ledd = 4. * np.pi * 1.67262158e-24 * c / 6.65245893699e-25 * GM
        eta = 1 - np.sqrt(8 / 9)
        Mcrit = Ledd / (c**2 * eta)
        fr = Freddi(wind=b'SS73C', windparams={}, Mx=Mx, Mdot0=Mcrit, Mdotout=Mcrit*Rout/Rin, initialcond=b'sinusF',
                    time=100*DAY, tau=1*DAY, rout=Rout)
        for state in fr:
            pass
        h = state.h
        F = Mcrit / Rin / (3 * GM) * (h**3 - h[0]**3)
        assert_allclose(state.F, F, rtol=1e-5)


class StationaryWindATestCase(unittest.TestCase):
    _k_A0 = 10.0

    @unittest.skipIf(scipy is None, 'No scipy module is available')
    def test(self):
        """Stationary disk with outflow rate proportional to accretion rate

        Outflow rate $d\dot M / dA ~ dF/dh$:

        $d\dot M_w / dA = -k_A (h - h_{in}) / (h_{out} - h_{in})^2
            \times 4\pi h^3 / (GM)^2 \times dF/dh$

        Stationary solution is
        $F ~ erfi(\sqrt{k_A/2} (h - h_{in})/(h_{out} - h_{in}))$

        The initial condition is non-stationary and does not satisfies outer
        boundary condition.

        """
        Mdot = 1e18
        fr = Freddi(wind=b'__testA__', windparams={'kA': self._k_A0}, Mdotout=Mdot, initialcond=b'sinusF', Mdot0=Mdot,
            time=1000*DAY, tau=1*DAY, Nx=10000, gridscale=b'linear')
        for state in fr:
            pass
        h = state.h
        F = (Mdot * (h[-1] - h[0]) * np.sqrt(np.pi/2/self._k_A0) * np.exp(-self._k_A0/2)
             * scipy.special.erfi(np.sqrt(self._k_A0/2) * (h - h[0]) / (h[-1] - h[0])))
        assert_allclose(state.F, F, rtol=1e-6)


class StationaryWindBTestCase(unittest.TestCase):
    _k_B0 = 16.0

    def test(self):
        """Stationary disk with outflow rate proportional to viscous torque

        Outflow rate $d\dot M / dA ~ F$:

        $d\dot M_w / dA = -k_B / (h_{out} - h_{in})^2 \times 4\pi h^3 / (GM)^2
            \times F$

        Stationary solution is
        $F ~ sinh(\sqrt{k_B} (h - h_{in})/(h_{out} - h_{in}))$

        The initial condition is non-stationary and does not satisfies outer
        boundary condition.

        """
        Mdot = 1e18
        fr = Freddi(wind=b'__testB__', windparams={'kB': self._k_B0}, Mdotout=Mdot, initialcond=b'sinusF', Mdot0=Mdot,
                    time=1000*DAY, tau=1*DAY, Nx=10000, gridscale=b'linear')
        for state in fr:
            pass
        h = state.h
        F = (Mdot / np.sqrt(self._k_B0) * (h[-1] - h[0])
             / np.sinh(np.sqrt(self._k_B0))
             * np.sinh((h - h[0]) / (h[-1] - h[0]) * np.sqrt(self._k_B0)))
        assert_allclose(state.F, F, rtol=1e-3)


class StationaryWindCTestCase(unittest.TestCase):
    _k_hw = 0.5
    _k_C0 = 3.0

    def test(self):
        """Stationary disk with constant outflow rate

        Outflow rate $d\dot M / dA$ doesn't depend on F nor dF/dh

        $d\dot M_w / dA = k_C / 2 (\cos(2\pi (h - h_w)/(h_{out} - h_w)) - 1)$
        for $h > h_w$ and zero otherwise

        The initial condition is non-stationary and does not satisfies outer
        boundary condition.

        """
        Mdot = 1e18
        fr = Freddi(wind=b'__testC__', windparams={'kC': self._k_C0}, Mdotout=Mdot, initialcond=b'sinusF', Mdot0=Mdot,
                    time=10000*DAY, tau=10*DAY, gridscale=b'linear')
        for state in fr:
            pass
        h = state.h
        h_w = self._k_hw * h[-1]
        C0 = self._k_C0 * Mdot / (h[-1] - h[0])
        F = (Mdot - C0 / 2 * (h[-1] - h_w)) * (h - h[0])
        i = h >= h_w
        F[i] += C0 / 2 * (
            -((h[-1] - h_w) / (2 * np.pi))**2 * (1 - np.cos(2*np.pi * (h[i] - h_w) / (h[-1] - h_w)))
            + (h[i] - h_w)**2 / 2
        )
        assert_allclose(state.F, F, rtol=1e-5)


class LipunovaShakuraTestCase(unittest.TestCase):
    """Test Lipunova-Shakura 2000 quasi-stationary solutions

    The analytical stationary solution is $F ~ (h - h_{in}) / h * P(h)$, where
    $P(h)$ is trinomial.

    The initial condition is non-stationary and does not satisfies outer
    boundary condition.

    """

    @staticmethod
    def F_LS(h, a, kl):
        xi = h / h[-1]
        f = a[0] * xi + a[1] * xi**kl[0] + a[2] * xi**kl[1]
        return f * (h - h[0]) / (h[-1] - h[0]) / xi

    def check_freddi(self, opacity, a, kl):
        fr = Freddi(opacity=opacity.encode(), time=1000*DAY, tau=1*DAY)
        for state in fr:
            pass
        h = state.h
        F = state.F
        assert_allclose(F / F[-1], self.F_LS(h, a, kl), rtol=1e-3)

    def test_kramers(self):
        """Lipunova-Shakura 2000 quasi-stationary solution for Kramers opacity"""
        opacity = 'Kramers'
        a = (1.430, -0.460, 0.030)
        kl = (3.5, 6.0)
        self.check_freddi(opacity, a, kl)

    def test_opal(self):
        """Lipunova-Shakura 2000 quasi-stationary solution for "Opal" opacity"""
        opacity = 'OPAL'
        a = (1.400, -0.425, 0.025)
        kl = (11. / 3., 19. / 3.)
        self.check_freddi(opacity, a, kl)


class ShieldPowerLawWindTestCase(unittest.TestCase):
    """Shield et al. 1986, IVb), eq. 4.2"""

    _k_C = 1
    _r_wind = 0.9

    def test_stationary_OPAL_q0(self):
        Mx = 2e34
        GM = 6.673e-8 * Mx
        Mdotout = 1e18
        rout = 1e11
        Mdotin = Mdotout / (1 + self._k_C)
        fr = Freddi(wind=b'__testC_q0_Shields1986__', windparams={'kC': self._k_C, 'Rwind': self._r_wind},
                    F0=Mdotin*np.sqrt(GM*rout), Mdotout=Mdotout, rout=rout, Mx=Mx,
                    initialcond=b'powerF', powerorder=1,
                    time=1000*DAY, tau=1*DAY, Nx=10000, gridscale=b'linear')
        for state in fr:
            pass
        h = state.h
        h_w = h[-1] * np.sqrt(self._r_wind)
        F = Mdotin * (h - h[0])
        i = h > h_w
        F[i] += self._k_C * Mdotin * (h[i] * np.log(h[i] / h_w) - (h[i] - h_w)) / np.log(h[-1] / h_w)
        assert_allclose(state.F, F, rtol=1e-3)
