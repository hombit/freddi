#!/usr/bin/env python3


import argparse
import re
import os, os.path
import numpy as np
from io import StringIO
from multiprocessing import Pool
from functools import partial

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




def kerrMdot(obs_filename=None, cloptions={}, start_from_peak=False):
    if obs_filename is None:
        obs_filename = '/Users/hombit/Dropbox/X-ray_novae_modeling (2) (1)/data_and_plots/Mdot-t/Min_simpl_kerrbb_gauss_smedge_ak0.9.v2.dat'
        Mx = 9.4
        kerr = 0.9
        spectrum_fit = 'kerrbb_gauss_smedge'
    else:
        basename = os.path.basename(obs_filename)
        match_result = re.match( r'Min_(?P<spectrum_fit>.+)_ak_(?P<kerr>0.\d+)_m(?P<Mx>\d+(\.\d+)?).dat', basename )
        match_dict = match_result.groupdict()
        Mx = float( match_dict['Mx'] )
        kerr = float( match_dict['kerr'] )
        spectrum_fit = match_dict['spectrum_fit']

    cloptions.update({
        'kerr' : kerr,   
    })

    Time_shift = 446
    Time_peak = 447

    obscolumns = ['Day', 'Mdot', 'Mdot_negerr', 'Mdot_poserr', 'tstart', 'tstop', 'chi2']
    obs = np.genfromtxt(obs_filename, names=obscolumns)
    if start_from_peak:
        obs = obs[ obs['Day'] > Time_peak ]
    obs = obs[ obs['chi2'] < 2 ] # Take dots with good statistics
    obs = obs[ obs['Mdot'] > 5e-3 ] # Do not take dots with too small Mdot
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
        mulF0_min=0.01,
        mulF0_max=2.0,
    )

    fred = FRED(
        obs=obs,
        sp=sp,
        op=op,
        absolute=True,
        flux_obs='Mdot',
        flux_model='Mdot',
        flux_model_func=lambda Mdot: Mdot / 1e18,
        flux_fit_model='None',
        cloptions=cloptions,
    )

    #alpha = 0.64902
    #fred.print_params( fred.fit_F0(alpha), alpha )

    # fred.print_params( *fred.fit_F0alpha() )

    with StringIO() as stream:
        #try:
            fred.print_params(
                *fred.fit_F0alpha(),
                #1.7642e+36, 0.53963,
                stream=stream,
                oneline=True,
                additional_fields={
                    'spectrum_fit' : spectrum_fit,
                }
            )
            line = stream.getvalue()
        #except Exception as e:
        #    line = '# ERROR: Mx = {}\tkerr = {}. {}\n'.format(Mx, kerr, e)
        #    print("Cannot compute something: "+line)

    return line


def process_kerrMdot(filenames, cloptions={}, start_from_peak=False, multiproc=True):
    partfunc = partial(kerrMdot, cloptions=cloptions, start_from_peak=start_from_peak)

    if multiproc and len(filenames) > 1:
        with Pool(1) as p:
            lines = p.map(partfunc, filenames)
    else:
        lines = map(partfunc, filenames)
 
    return lines


#############



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find your value of alpha in 1543!')
    parser.add_argument(
        'filename',
        action = 'store',
        nargs = '+',
        help = 'path to Mdot-t file'
    )
    parser.add_argument(
        '-C', '--Cirr',
        dest = 'Cirr',
        action = 'store',
        default = 1e-3,
        type = float,
        help = 'parameter of irradiation',
    )
    parser.add_argument(
        '--results',
        dest = 'results',
        action = 'store',
        default = 'results.dat',
        help = 'where to output results of optimization',
    )
    args = parser.parse_args()
    Cirr = args.Cirr
    filenames = args.filename

    cloptions = {
        'initialcond' : 'quasistat',
#        'powerorder' : 6,
        'Thot' : 1e4,
        'boundcond' : 'Tirr',
        'Cirr' : Cirr,
        'distance' : 4.937,
        'Nx' : 1000,
        'gridscale' : 'linear',
        'opacity' : 'OPAL',
    }

    lines = process_kerrMdot(filenames, cloptions, start_from_peak=True, multiproc=True)
    with open(args.results, 'w') as f:
        for line in lines:
            f.write(line)
