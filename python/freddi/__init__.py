from functools import partial, singledispatch

import numpy as np

from ._freddi import _Freddi, _FreddiNeutronStar
from .evolution_result import EvolutionResult


class _MetaFreddi(type(_Freddi)):
    def __new__(mcs, name, bases, attrs, boost_cls=object):
        attrs['_boost_class'] = boost_cls
        test_obj = boost_cls(**boost_cls._required_args())
        for attr_name in dir(boost_cls):
            if attr_name.startswith('_'):
                continue
            cls_attr = getattr(boost_cls, attr_name)
            if not isinstance(cls_attr, property):
                continue
            attrs[attr_name] = property(partial(mcs._boost_cls_property, attr_name))
            obj_attr = getattr(test_obj, attr_name)
            if isinstance(obj_attr, np.ndarray):
                attrs['first_' + attr_name] = property(partial(mcs._first_last_getter, 'first', attr_name))
                attrs['last_' + attr_name] = property(partial(mcs._first_last_getter, 'last', attr_name))

        return super().__new__(mcs, name, bases, attrs)

    def __init__(cls, name, bases, attrs, boost_cls=object):
        super().__init__(name, bases, attrs)

    @staticmethod
    def _boost_cls_property(attr_name, self):
        return getattr(self._freddi, attr_name)

    @staticmethod
    def _first_last_getter(first_or_last, attr_name, self):
        return getattr(self._freddi, attr_name)[getattr(self._freddi, first_or_last)]


class _BasePyFreddi:
    def __init__(self, **kwargs):
        self._freddi = self._boost_class(**kwargs)
        self._kwargs = self._freddi._kwargs

    @classmethod
    def _from_boost(cls, boost_freddi):
        obj = cls.__new__(cls)
        obj._freddi = boost_freddi
        return obj

    @classmethod
    def from_astropy(cls, **kwargs):
        """Alternative Freddi constructor accepting astropy Quantity

        All physical values can be passed as `astropy.units.Quantity`, and all
        string values can be passed as `str`. Requires `astropy` package

        Parameters
        ----------
        **kwargs :
            `Freddi` arguments

        Returns
        -------
        Freddi

        """
        from collections.abc import Mapping
        from functools import singledispatch
        
        from astropy.units import Quantity

        @singledispatch
        def convert(value):
            return value

        @convert.register(Quantity)
        def _(value):
            return value.cgs.value

        @convert.register(Mapping)
        def _(value):
            return {k: convert(v) for k, v in value.items()}

        kwargs = {key: convert(value) for key, value in kwargs.items()}
        
        return cls(**kwargs)

    @classmethod
    def alt(cls, **kwargs):
        """Alias to .from_astropy()"""
        return cls.from_astropy(**kwargs)

    def evolve(self):
        """Calculate disk evolution

        Returns
        -------
        EvolutionResults

        """
        return EvolutionResult(self)

    def _flux_hot(self, lmbd, phase):
        del phase
        return self._freddi._flux_hot(lmbd)

    def _flux_cold(self, lmbd, phase):
        del phase
        return self._freddi._flux_cold(lmbd)

    def _flux_star(self, lmbd, phase):
        if phase is None:
            raise ValueError('Phase must be specified if star flux is required')
        return self._freddi._flux_star(lmbd, phase)

    def flux(self, lmbd, region='hot', phase=None):
        region = region.lower()
        if 'hot'.startswith(region):
            flux = self._flux_hot
        elif 'cold'.startswith(region):
            flux = self._flux_cold
        elif 'disk'.startswith(region):
            def flux(lmbd, phase):
                return self._flux_hot(lmbd, phase) + self._flux_cold(lmbd, phase)
        elif 'star'.startswith(region):
            flux = self._flux_star
        elif 'all'.startswith(region):
            def flux(lmbd, phase):
                return self._flux_hot(lmbd, phase) + self._flux_cold(lmbd, phase) + self._flux_star(lmbd, phase)
        else:
            raise ValueError(f'Zone {region} is not supported')

        lmbd = np.asarray(lmbd)
        arr = np.empty_like(lmbd, dtype=float)
        for i, x in np.ndenumerate(lmbd):
            arr[i] = flux(float(x), phase)
        return arr

    def __iter__(self):
        for value in self._freddi:
            value = self._from_boost(value)
            yield value


class Freddi(_BasePyFreddi, metaclass=_MetaFreddi, boost_cls=_Freddi):
    pass


class FreddiNeutronStar(_BasePyFreddi, metaclass=_MetaFreddi, boost_cls=_FreddiNeutronStar):
    pass


__all__ = ('Freddi', 'FreddiNeutronStar')
