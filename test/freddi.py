import unittest

from freddi import Freddi


class ChangeArgsTestCase(unittest.TestCase):
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
