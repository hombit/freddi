import unittest

from numpy.testing import assert_allclose, assert_array_equal

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


class NSMdotFractionTestCase(unittest.TestCase):
    """fp(r/r_cor) function"""
    basic_kwargs = dict(
        freqx=500,  # Hz
        rout=1e11,  # cm
        Bx=1e8,  # G
        Rdead=1e7,  # cm
        time=100 * 86400,  # s
        tau=0.1 * 86400,  # s
    )

    def test_no_outflow(self):
        fr = FreddiNeutronStar(fptype='no-outflow', **self.basic_kwargs)
        assert_array_equal(1, fr.fp)
