import numpy as np


class EvolutionResult:
    """Temporal distribution of various disk parameters

    Parameters
    ----------
    states: Iterable of states, e.g. list of Freddi
        States to be stored as results

    Methods
    -------
    flux(lmbd, region, phase) : array
        Optical flux of the disk

    """
    def __init__(self, states):
        self.__states = tuple(states)
        self.__len_states = len(self.__states)

    def flux(self, lmbd, region='hot', phase=None) -> np.ndarray:
        """Optical flux of the disk

        Parameters
        ----------
        lmbd : array_like
            Observation wavelength
        region: str, optional
            Object optically luminous component, could be one of: `hot`, `cold`,
            `disk`, `star`, `all`
        phase: float or None, optional
            Phase of the observation, must be specified if `region` is `star` or
            `all`

        Returns
        -------
        array

        """
        lmbd = np.asarray(lmbd)
        arr = np.empty((self.__len_states,) + lmbd.shape, dtype=np.float)
        for i in range(self.__len_states):
            arr[i] = self.__states[i].flux(lmbd, region, phase)
        return arr

    def __getattr__(self, attr) -> np.ndarray:
        first_val = np.asarray(getattr(self.__states[0], attr))
        if first_val.ndim == 0:
            arr = np.empty(self.__len_states, dtype=first_val.dtype)
            arr[0] = first_val
            for i in range(1, self.__len_states):
                arr[i] = getattr(self.__states[i], attr)
            return arr
        elif first_val.ndim == 1:
            arr = np.full((self.__len_states, self.__states[0].Nx), np.nan, dtype=first_val.dtype)
            for i in range(self.__len_states):
                value = np.asarray(getattr(self.__states[i], attr))
                arr[i, self.__states[i].first:self.__states[i].last+1] = value
            return arr
        else:
            raise ValueError("{}.ndim > 1 ({})".format(attr, first_val.ndim))
