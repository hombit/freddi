#!/usr/bin/env python3


import numpy as np
from functools import reduce


def lcmean( source, dest, interval=1 ):
    names = ['TJD', 'DoY', 'DaP', 'flux', 'err']
    s = np.genfromtxt( source, skip_header=10, names=names, usecols=tuple(range(5)) )
    
    tjd = s['TJD'][0]
    interval_set = []
    d = np.zeros( s.shape, dtype=s.dtype )
    n_d = 0
    for x in s:
        if x['TJD'] < tjd + interval:
            interval_set.append(x)
        else:
            n_set = len(interval_set)
            if n_set == 1:
                d[n_d] = x
            else:
                interval_set = np.array( [interval_set], dtype=s.dtype )
                for name in names[:3]:
                    d[name][n_d] = interval_set[name].mean()
                sigma2 = interval_set['err']**2
                p = ( 1. / sigma2 ).sum()
                d['flux'][n_d] = ( interval_set['flux'] / sigma2 ).sum() / p
                d['err'][n_d] = sqrt(1. / p) + np.sqrt( ( (interval_set['flux'] - d['flux'][n_d])**2 / sigma2 ).sum() / p / (n_set-1) )
            n_d += 1
            tjd = np.ceil( x['TJD'] )
            interval_set = [x]
            
    with open(dest, 'w') as fd:
        fd.write(
            '# ' + reduce( lambda x,y: '{} {}'.format(x,y), names ) + '\n'
        )
        for i in range(n_d):
            fd.write(
                reduce( lambda x,y: '{} {}'.format(x,y), d[i] ) + '\n'
            )
