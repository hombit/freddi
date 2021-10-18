from configparser import ConfigParser
import os

from freddi import Freddi, FreddiNeutronStar


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

DAY = 86400

# Should contain all dimensional qu
_UNITS = dict(
    Mx=1.98892e33,
    Mopt=1.98892e33,
    period=DAY,
    distance=1000 * 3.08567758135e18,
    time=DAY,
    tau=DAY,
)


def _to_float(x):
    try:
        return float(x)
    except ValueError:
        return x


def default_freddi_kwargs():
    header = 'DEFAULT'
    with open(os.path.join(DATA_DIR, 'freddi.ini')) as f:
        freddi_ini = '[{}]\n'.format(header) + f.read()

    cp = ConfigParser(inline_comment_prefixes=['#'])
    cp.optionxform = str
    cp.read_string(freddi_ini)
    kwargs = dict(cp[header].items())
    kwargs = {name: _to_float(value) for name, value in kwargs.items()}
    kwargs.update({name: unit * kwargs[name] for name, unit in _UNITS.items() if name in kwargs})
    return kwargs


def defaulted_kwargs(**kwargs):
    return {**default_freddi_kwargs(), **kwargs}


def freddi_w_default(**kwargs):
    return Freddi(**defaulted_kwargs(**kwargs))


def freddi_ns_w_default(**kwargs):
    return FreddiNeutronStar(**defaulted_kwargs(**kwargs))
