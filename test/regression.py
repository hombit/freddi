import glob
import os
import re
import unittest
from configparser import ConfigParser

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from freddi import Arguments, Freddi



DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class RegressionTestCase(unittest.TestCase):
    def setUp(self):
        self.data_files = glob.glob(os.path.join(DATA_DIR, '*.dat'))

    @staticmethod
    def load_config(data_file, arguments_to_remove=()):
        section = 'default'
        section_str = '[{}]\n'.format(section)
        regexp = re.compile(r'^# ([^-])')
        with open(data_file) as f:
            lines = (regexp.sub(r'\1', line) for line in f if regexp.match(line))
            config_str = section_str + ''.join(lines)
        config = ConfigParser(default_section=None)
        config.optionxform = str
        config.read_string(config_str)
        config = dict(config[section])
        for kwarg in arguments_to_remove:
            config.pop(kwarg, None)
        for name, value in config.items():
            try:
                config[name] = float(value)
            except ValueError:
                config[name] = value.encode()
        config['cgs'] = False
        return config

    def test(self):
        columns_to_compare = ('Mdot', 'Mdisk', 'Lx', 'mU', 'mB', 'mV', 'mR', 'mI', 'mJ')
        for data_file in self.data_files:
            config = self.load_config(data_file, ('dir', 'prefix', 'fulldata'))
            args = Arguments(**config)
            f = Freddi(args)
            result = f.evolve()
            data = np.genfromtxt(data_file, names=True)
            assert_equal(result.i_t[:, 0], np.arange(f.Nt + 1))
            assert_equal(result.t[:, 0], data['t'] * 86400)
            for column in columns_to_compare:
                assert_allclose(getattr(result, column)[:, 0], data[column], rtol=1e-5,
                                err_msg='Column {}:'.format(column))
