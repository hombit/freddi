import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from freddi import FreddiNeutronStar

from .util import default_freddi_kwargs

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
    freddi_kwargs = {**default_freddi_kwargs(), 'inversebeta': inverse_beta, 'Bx': Bx, 'freqx': freqx, 'Rx': Rx}

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
    Rx = 1e6  # cm
    rin = 1.3e6  # cm larger Rx and Risco
    basic_kwargs = {
        **default_freddi_kwargs(),
        'Mdot0':1e16,
        'initialcond':'quasistat',
        'Mx':1.4 * 2e33,  # g
        'Rx':Rx,  # cm
        'freqx':500,  # Hz
        'rin': rin,  # cm
        'rout':1e11,  # cm
        'Bx':1e8,  # G
        'Rdead':5e6,  # cm
        'time':200 * 86400,  # s
        'tau':0.1 * 86400,  # s
    }

    def test_no_outflow(self):
        fr = FreddiNeutronStar(fptype='no-outflow', **self.basic_kwargs)
        result = fr.evolve()
        assert_array_equal(1, result.fp)

    def test_propeller(self):
        fr = FreddiNeutronStar(fptype='propeller', **self.basic_kwargs)
        result = fr.evolve()
        assert_array_equal(0, result.fp)

    def test_propeller_unity_fp(self):
        """fp must be unity if R_in <= Rx or R_in <= Risco"""
        kwargs = self.basic_kwargs.copy()
        kwargs['rin'] = self.Rx
        fr = FreddiNeutronStar(fptype='propeller', **kwargs)
        result = fr.evolve()
        idx = fr.first_R == self.Rx
        assert_array_equal(1, result.fp[idx])

    def test_corotation_block(self):
        fr = FreddiNeutronStar(fptype='corotation-block', **self.basic_kwargs)
        result = fr.evolve()
        r = result.first_R / fr.R_cor
        assert np.any(r < 1)
        assert np.any(r > 1)
        assert_array_equal(r < 1, result.fp)

    def test_geometrical_chi_zero(self):
        fr_corotation_block = FreddiNeutronStar(fptype='corotation-block', **self.basic_kwargs)
        result_corotation_block = fr_corotation_block.evolve()
        fr_geometrical = FreddiNeutronStar(fptype='geometrical', fpparams={'chi': 0}, **self.basic_kwargs)
        result_geometrical = fr_geometrical.evolve()
        assert_array_equal(result_corotation_block.fp, result_geometrical.fp)

    def test_geometrical_chi_nonzero(self):
        fr = FreddiNeutronStar(fptype='geometrical', fpparams={'chi': 30}, **self.basic_kwargs)
        result = fr.evolve()
        fp = result.fp
        r = result.first_R / fr.R_cor
        assert np.all(fp >= 0) and np.all(fp <= 1)
        assert r[0] < 1 and fp[0] == 1, 'Initial inner radius should be less than corotation one and fp should equal zero'
        assert np.all(fp[r > 1] < 1)
