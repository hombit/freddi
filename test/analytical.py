import unittest

import numpy as np
from numpy.testing import assert_allclose

from freddi import Freddi


DAY = 86400


class ShakuraSunyaevSubctriticalTestCase(unittest.TestCase):
    def test(self):
        Mdot = 1e18
        fr = Freddi(initialcond=b'sinusF', Mdot0=Mdot/2, Mdotout=Mdot, time=1000*DAY, tau=1*DAY)
        for state in fr:
            pass
        h = state.h
        assert_allclose(state.F, Mdot * (h - h[0]))

class StationaryWindCTestCase(unittest.TestCase):
    _k_hw = 0.5
    _k_C0 = 3.0

    def test(self):
        Mdot = 1e18
        fr = Freddi(wind=b'__testC__', Mdotout=Mdot, initialcond=b'sinusF', Mdot0=Mdot,
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
    def F_LS(self, h, a, kl):
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
        opacity = 'Kramers'
        a = (1.430, -0.460, 0.030)
        kl = (3.5, 6.0)
        self.check_freddi(opacity, a, kl)

    def test_opal(self):
        opacity = 'OPAL'
        a = (1.400, -0.425, 0.025)
        kl = (11. / 3., 19. / 3.)
        self.check_freddi(opacity, a, kl)
