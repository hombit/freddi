import glob
import os
import re
import unittest
from collections import UserDict
from configparser import ConfigParser, RawConfigParser

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from freddi import Arguments, Freddi


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class RegressionTestCase(unittest.TestCase):
    def setUp(self):
        self.data_files = glob.glob(os.path.join(DATA_DIR, '*.dat'))

    class MultiDict(UserDict):
        def __setitem__(self, key, value):
            if isinstance(value, list):
                value = value[0]
            if key not in self.data:
                self.data[key] = value
                return
            if value == self.data[key]:
                return
            if isinstance(self.data[key], tuple):
                if value in self.data[key]:
                    return
                self.data[key] += (value,)
                return
            self.data[key] = (self.data[key], value)

    def convert_config_value(self, value):
        if isinstance(value, tuple):
            return [self.convert_config_value(v) for v in value]
        try:
            return float(value)
        except ValueError:
            return value.encode()

    def load_config(self, data_file, arguments_to_remove=()):
        section = 'section'
        section_str = '[{}]\n'.format(section)
        regexp = re.compile(r'^# ([^-])')
        with open(data_file) as f:
            lines = (regexp.sub(r'\1', line) for line in f if regexp.match(line))
            config_str = section_str + ''.join(lines)
        config = RawConfigParser(strict=False, dict_type=self.MultiDict, inline_comment_prefixes=('#', ';'))
        config.optionxform = str
        config.read_string(config_str)
        config = dict(config[section].items())
        for kwarg in arguments_to_remove:
            config.pop(kwarg, None)
        for name, value in config.items():
            config[name] = self.convert_config_value(value)
        if 'lambda' in config:
            config['lambdas'] = config.pop('lambda')
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
                                err_msg='File: {}, column {}:'.format(data_file, column))
            for i_lmbd, lmbd in enumerate(args.lambdas):
                column = 'Fnu{}'.format(i_lmbd)
                if column in data.dtype.names:
                    assert_allclose(result.flux(lmbd)[:, 0], data[column], rtol=1e-5,
                                    err_msg='File: {}, lambda {}, column {}:'.format(data_file, lmbd, column))
                else:
                    break
