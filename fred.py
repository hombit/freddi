#!/usr/bin/env python3


import re
import os, os.path
from glob import glob
from io import StringIO
from multiprocessing import Pool
import numpy as np

import fredlib
from fredlib import FRED, Tmax_SS, Tin_color_diskbb


def parkTin():
    obs_filename = 'obs/park1543.dat'
    obscolumns = ['ObsN', 'Date', 'MJD', 'ObsTime', 'Tin', 'Tin_poserr', 'Tin_negerr', 'PhotonIndex', 'PhotonIndex_poserr', 'PhotonIndex_negerr', 'PowerLawNorm', 'PowerLawNorm_poserr', 'PowerLawNorm_negerr', 'FeEdge', 'FeEdge_poserr', 'FeEdge_negerr', 'tauFe', 'tauFe_poserr', 'tauFe_negerr', 'Lx', 'Lx_poserr', 'Lx_negerr', 'chi2nu']
    obs = np.genfromtxt(obs_filename, names=obscolumns)
    obs = fredlib.rec_append_fields( obs, 'DaP', obs['MJD'] - 52445.5 )
    obs = fredlib.rec_append_fields( obs, 'err', 0.5 * (obs['Tin_poserr'] + obs['Tin_negerr']) )

    sp = fredlib.SystemParameters(
        Mx=6,
        Mopt=2.5,
        Period=1.1164,
        r_out=1,
        fcol=1.7,
    )

    op = fredlib.OptimizeParameters(
        Time=20,
        tau=0.2,
        t0_range=3,
        alpha_min= 0.1,
        alpha_max= 1.0,
        alpha_step=0.025,
        mulF0_min=0.1,
        mulF0_max=15.0,
    )

    fred = FRED(
        obs=obs,
        sp=sp,
        op=op,
        absolute=True,
        flux_obs='Tin',
        flux_model='Mdot',
        flux_model_func=lambda Mdot: Tin_color_diskbb( Tmax_SS(sp.Mx, Mdot, 6.*sp.Rg, units='keV'), fcol=sp.fcol ),
        flux_fit_model='t0',
        cloptions=(
            '--initialcond=power',
            '--powerorder=6',
            '--Thot=-1',
            '--Cirr=0.00',
    #        '--fulldata'
        ),
    )
    #fred.print_params( 1.9764e+38, 0.71178 )

    #fred.fit_alpha()

    #alpha = 0.71178
    #fred.print_params( fred.fit_F0(alpha), alpha )

    fred.print_params( 1.4893e+37, 0.40401 )
    print( fred.find_errors( 1.4893e+37, 0.40401, sigma_level=3 ) )




def kerrMdot(obs_filename=None):
    if obs_filename is None:
        obs_filename = '/Users/hombit/Dropbox/X-ray_novae_modeling (2) (1)/data_and_plots/Mdot-t/Min_simpl_kerrbb_gauss_smedge_ak0.9.v2.dat'
        Mx = 9.4
        kerr = 0.9
        spectrum_fit = 'kerrbb_gauss_smedge'
    else:
        basename = os.path.basename(obs_filename)
        match_result = re.match( r'Min_(?P<spectrum_fit>.+)_ak_(?P<kerr>0.\d+)_m(?P<Mx>\d+).dat', basename )
        match_dict = match_result.groupdict()
        Mx = float( match_dict['Mx'] )
        kerr = float( match_dict['kerr'] )
        spectrum_fit = match_dict['spectrum_fit']
    
    Time_shift = 445

    obscolumns = ['Day', 'Mdot', 'Mdot_negerr', 'Mdot_poserr', 'tstart', 'tstop', 'chi2']
    obs = np.genfromtxt(obs_filename, names=obscolumns)
    obs = obs[ obs['chi2'] < 2 ]
    obs = fredlib.rec_append_fields( obs, 'DaP', obs['Day'] - Time_shift )
    obs = fredlib.rec_append_fields( obs, 'err', 0.5 * (np.abs(obs['Mdot_poserr']) + np.abs(obs['Mdot_negerr'])) )

    sp = fredlib.SystemParameters(
        Mx=Mx,
        Mopt=2.5,
        Period=1.1164,
        r_out=None,
        fcol=1.7,
    )

    op = fredlib.OptimizeParameters(
        Time=35,
        Time_shift=Time_shift,
        tau=0.25,
        t0_range=3,
        alpha_min= 0.1,
        alpha_max= 2.0,
        alpha_step=0.025,
        mulF0_min=1.0,
        mulF0_max=15.0,
    )

    fred = FRED(
        obs=obs,
        sp=sp,
        op=op,
        absolute=True,
        flux_obs='Mdot',
        flux_model='Mdot',
        flux_model_func=lambda Mdot: Mdot / 1e18,
        flux_fit_model='t0',
        cloptions={
            'initialcond' : 'power',
            'powerorder' : 6,
            'Thot' : 1e4,
#            'boundSigma' : None,
            'Cirr' : 1e-4,
            'kerr' : kerr,
            'distance' : 4.937,
            'Nx' : 1000,
            'gridscale' : 'linear',
            'opacity' : 'OPAL',
        },
    )

    #alpha = 0.64902
    #fred.print_params( fred.fit_F0(alpha), alpha )

    # fred.print_params( *fred.fit_F0alpha() )

    with StringIO() as stream:
        try:
            fred.print_params(
                *fred.fit_F0alpha(),
                #2.0522e38, 0.76485,
                stream=stream,
                oneline=True,
                additional_fields={
                    'spectrum_fit' : spectrum_fit,
                }
            )
            line = stream.getvalue()
        except Exception as e:
            line = '# ERROR: Mx = {}\tkerr = {}. {}\n'.format(Mx, kerr, e)
            print("Cannot compute something: "+line)

    return line


def process_kerrMdot(filenames, multiproc=True):
    if multiproc:
        with Pool() as p:
            lines = p.map( kerrMdot, filenames )
    else:
        lines = [ line for line in map(kerrMdot, filenames) ]
    
    return lines


#############



if __name__ == '__main__':
    #print(kerrMdot( '/Users/hombit/Dropbox/X-ray_novae_modeling (2) (1)/data_and_plots/Mdot-t/Min_simpl_kerrbb_laor_smedge_ak_0.0_m6.dat' ))
    
    # parkTin()
    
    from sys import argv
    if len(argv) < 2:
        filenames = ['/Users/hombit/Dropbox/X-ray_novae_modeling (2) (1)/data_and_plots/Mdot-t/Min_simpl_kerrbb_laor_smedge_ak_0.0_m6.dat']
    else:
        filenames = argv[1:]

    if len(filenames) == 1:
        print( kerrMdot(filenames[0]) )
    else:
        lines = process_kerrMdot(filenames, multiproc=True)
        with open('results.dat', 'w') as f:
            for line in lines:
                f.write(line)
