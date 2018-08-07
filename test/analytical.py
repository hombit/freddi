import unittest

from numpy.testing import assert_allclose

from freddi import Freddi


class StationaryDiskTestCase(unittest.TestCase):
    def test_shakura_sunayev_subcritical(self):
        Mdot = 1e18
        day = 86400
        fr = Freddi(initialcond=b'sinusF', Mdot0=Mdot/2, Mdotout=Mdot, time=1000*day, tau=1*day)
        for state in fr:
            pass
        h = state.h
        assert_allclose(state.F, Mdot * (h - h[0]))
