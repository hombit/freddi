import unittest

from numpy.testing import assert_allclose

from freddi import Freddi


DAY = 86400


class StationaryDiskTestCase(unittest.TestCase):
    def test_shakura_sunayev_subcritical(self):
        Mdot = 1e18
        fr = Freddi(initialcond=b'sinusF', Mdot0=Mdot/2, Mdotout=Mdot, time=1000*DAY, tau=1*DAY)
        for state in fr:
            pass
        h = state.h
        assert_allclose(state.F, Mdot * (h - h[0]))

    def test_lipunova_shalura(self):
        fr = Freddi(opacity=b'Kramers', time=1000*DAY, tau=1*DAY)
        for state in fr:
            pass
        h = state.h
        F = state.F
        xi = h / h[-1]
        a0 = 1.430
        a1 = -0.460
        a2 = 0.030
        k = 3.5
        l = 6.0
        f = a0 * xi + a1 * xi**k + a2 * xi**l
        F_LS = F[-1] * f * (h - h[0]) / (h[-1] - h[0]) * h[-1] / h
        assert_allclose(F[-10:], F_LS[-10:], rtol=1e-3)
