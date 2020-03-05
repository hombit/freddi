import glob
import os
import re
import unittest
from collections import UserDict
from configparser import ConfigParser, RawConfigParser

import numpy as np
from numpy.testing import assert_allclose, assert_equal
from parameterized import parameterized

from freddi import Freddi

from .util import DATA_DIR


class RegressionTestCase(unittest.TestCase):
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
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value.encode()

    def load_config(self, data_file, arguments_to_remove=()):
        section = 'section'
        section_str = '[{}]\n'.format(section)
        with open(data_file) as f:
            content = f.read()
        parameters = re.search('### Parameters.+?###', content, flags=re.MULTILINE | re.DOTALL).group()
        parameters = re.sub('^# ', '', parameters, flags=re.MULTILINE)
        config_str = section_str + parameters
        config = RawConfigParser(strict=False, dict_type=self.MultiDict, inline_comment_prefixes=('#', ';'))
        config.optionxform = str
        config.read_string(config_str)
        config = dict(config[section].items())
        for kwarg in arguments_to_remove:
            config.pop(kwarg, None)
        for name, value in config.items():
            config[name] = self.convert_config_value(value)
        config['__cgs'] = False
        return config

    columns_to_compare = ('Mdot', 'Mdisk', 'Lx', 'mU', 'mB', 'mV', 'mR', 'mI', 'mJ')

    @parameterized.expand(glob.glob(os.path.join(DATA_DIR, '*.dat')))
    def test(self, data_file):
        config = self.load_config(data_file, ('dir', 'prefix', 'fulldata', 'precision', 'config'))
        lmbd_ = np.array(config.pop('lambda', [])) * 1e-8
        f = Freddi(**config)
        result = f.evolve()
        data = np.genfromtxt(data_file, names=True)
        with self.subTest('Number of lines'):
            assert_equal(result.i_t, np.arange(f.Nt + 1))
        with self.subTest('Time'):
            assert_equal(result.t, data['t'] * 86400)
        for column in self.columns_to_compare:
            with self.subTest(column):
                assert_allclose(getattr(result, column), data[column], rtol=1e-4,
                                err_msg='File: {}, column {}:'.format(data_file, column))
        for i_lmbd, lmbd in enumerate(lmbd_):
            column = 'Fnu{}'.format(i_lmbd)
            if column in data.dtype.names:
                with self.subTest(column):
                    assert_allclose(result.flux(lmbd), data[column], rtol=1e-4,
                                    err_msg='File: {}, lambda {}, column {}:'.format(data_file, lmbd, column))
            else:
                break
