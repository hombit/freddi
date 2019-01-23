import unittest

import numpy as np
from numpy.testing import assert_allclose

from freddi import FreddiNeutronStar

try:
    import scipy.interpolate
except ImportError:
    scipy = None


class FmagnDerivativesTestCase(unittest.TestCase):
    """Integrate first and second derivatives of magnetic torque"""
    inverse_beta = 0.5  # not zero and not unity
    Bx = 1e8  # G
    freqx = 100  # Hz
    Rx = 1e6  # cm
    freddi_kwargs = dict(inversebeta=inverse_beta, Bx=Bx, freqx=freqx, Rx=Rx)

    @unittest.skipIf(scipy is None, 'No scipy module is available')
    def test_first_derivative(self):
        """Integral of dFmagn_dh == Fmagn"""
        fr = FreddiNeutronStar(**self.freddi_kwargs)
        h = fr.h.copy()
        spline = scipy.interpolate.CubicSpline(h, fr.dFmagn_dh.copy())
        Fmagn = spline.antiderivative(1)(h)
        assert_allclose(fr.Fmagn - fr.Fmagn[0], Fmagn, rtol=1e-5, atol=0)

    @unittest.skipIf(scipy is None, 'No scipy module is available')
    def test_second_derivative(self):
        """Integral of dF2magn_dh2 == dFmagn_dh"""
        fr = FreddiNeutronStar(**self.freddi_kwargs)
        h = fr.h.copy()
        spline = scipy.interpolate.CubicSpline(h, fr.d2Fmagn_dh2.copy())
        dFmagn_dh = spline.antiderivative(1)(h)
        assert_allclose(fr.dFmagn_dh - fr.dFmagn_dh[0], dFmagn_dh, rtol=1e-5, atol=0)
