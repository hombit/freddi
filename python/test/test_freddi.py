import unittest

import numpy as np

from freddi import Freddi


class ChangeArgsTestCase(unittest.TestCase):
    @unittest.skip('Changing arguments is not supported any more')
    def test_Cirr(self):
        fr = Freddi(Cirr=1e-3)
        self.assertGreater(fr.last_Tirr, 0)
        fr.Cirr = 0
        self.assertEqual(fr.Cirr, 0)
        for state in fr:
            pass
        self.assertEqual(state.last_Cirr, 0)
        self.assertEqual(state.last_Tph, state.last_Tph_vis)
        self.assertEqual(state.last_Tirr, 0)


class TwoDimensionalValueTestCase(unittest.TestCase):
    def test(self):
        Nx = 1000
        fr = Freddi(Mx=1e34, Mopt=1e33, period=2e4,
                    F0=2e38, Thot=1e4, initialcond='sineF',
                    alpha=0.25, distance=1e19,
                    Nx=Nx, time=50 * 86400)
        evolution_result = fr.evolve()

        np.testing.assert_equal(evolution_result.Nx, Nx)

        nan_idx = np.vstack([
            (
                [True] * f
                + [False] * (l - f + 1)
                + [True] * (Nx - l - 1)
            ) for f, l in zip(evolution_result.first, evolution_result.last)
        ])
        for attr in ('R', 'h', 'F', 'Tph',):
            with self.subTest(attr):
                val = getattr(evolution_result, attr)
                self.assertTrue(np.all(np.isnan(val[nan_idx])))
                self.assertFalse(np.any(np.isnan(val[~nan_idx])))
