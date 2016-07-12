#!/usr/bin/env python3


import astropy.units
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing, subprocess
import os, os.path
import sys
from astropy import constants
from numpy.lib.recfunctions import rec_append_fields
from scipy.interpolate import interp1d
from scipy.optimize import brentq, curve_fit, differential_evolution
from scipy import special, stats


eV = astropy.units.eV.in_units(astropy.units.erg)
keV = 1000. * eV

DAY = 86400
k_B = constants.k_B.cgs.value
sigma_sb = constants.sigma_sb.cgs.value
G = constants.G.cgs.value
c = constants.c.cgs.value

Msun = constants.M_sun.cgs.value
Rsun = constants.R_sun.cgs.value

GMsun = G * Msun
Rgsun = GMsun / c**2

Mdot_edd = 1.7e18    


def Tmax_SS( Mx, Mdot, r_in, units='Kelvin' ):
    'Ref: Book "Accretion processes in astrophysics", chapter 1 (G. Lipunova)'

    T_kelvin = 2**0.75 * (3./7.)**1.75 * ( Mx*GMsun * Mdot / ( np.pi * sigma_sb * r_in**3. ) )**0.25
    if   units == 'Kelvin':
        return T_kelvin
    elif units == 'keV':
        return k_B * T_kelvin / keV
    else:
        raise ValueError( "Unknown units {}".format(units) )


def Tin_color_diskbb( TmaxSS, fcol=1.7 ):
    'Ref: Zimmerman et. al., AJ 618:832, 2005'
    return 1.05 * fcol * TmaxSS


def chi2(f, x, y, sigma=None):
    if sigma is None  or  (sigma == 0).any():
        sigma = np.ones_like(y)
    return (( (f(x) - y) / sigma )**2).sum() / (y.shape[0]-2)

def chi2_errors( chi2_func, x0, bounds, dof, sigma_level=1 ):
    chi2_quantile = stats.distributions.chi2.ppf( special.erf(sigma_level/np.sqrt(2)), df=dof ) / dof
    chi2_min = chi2_func(x0)
    if chi2_quantile > chi2_min:
        chi2_quantile = chi2_min
    print( chi2_min, chi2_quantile, sigma_level**2*chi2_min/chi2_quantile/dof)

    errors = np.empty( ( len(x0), 2 ) )
    for i in range( len(x0) ):
        bounds_trunc = list(bounds)
        del bounds_trunc[i]

        def chi2_func_trunc(x, x_i):
            x_list = list(x)
            x_list.insert(i, x_i)
            return chi2_func(x_list)

        search_func = lambda x_i: differential_evolution(
            lambda x: chi2_func_trunc(x, x_i) - chi2_min - sigma_level**2*chi2_min/chi2_quantile/dof,
            bounds_trunc,
            tol=2e-2,
            popsize=7,
            seed=np.random.mtrand.RandomState(),
        ).fun

        for j in range(2):
            errors[i,j] = abs(
                brentq(
                        search_func,
                        bounds[i][j], x0[i],
                        rtol=2e-2,
                    )
                - x0[i]
            )
    return errors

def disksize(Mx, Mopt, Period):
    Mass = (Mx+Mopt) * Msun
    Period_s = Period * DAY
    a = ( G * Mass * Period_s**2. / (2.*np.pi)**2. )**(1./3.) / Rsun
    q = (Mx / Mopt)**(1./3.)
    roche = (a * 0.49 * q**2.) / ( 0.6 * q**2. + np.log(1.+q) )
    return 0.8 * roche


class SystemParameters(object):
    def __init__(
        self,
        Mx,
        Mopt=None,
        Period=None,
        r_out=None,
        Mdot_max=None,
        fcol=1.7,
    ):
        self.Mx = Mx  
        self.Mopt = Mopt
        self.Period = Period
        self.Mdot_max = Mdot_max or Mdot_edd * Mx
        self.Rg = Rgsun * Mx
        self.fcol = fcol
        
        if not(Mopt is None) and not(Period is None) and (r_out is None):
            self.r_out = disksize( Mx, Mopt, Period )
        elif not(r_out is None):
            self.r_out = r_out
        else:
            raise RuntimeError( "You should use set at least a pair of Mopt and Period or r_out" )



class OptimizeParameters(object):
    def __init__(
        self,
        Time=50,
        Time_shift = 0.,
        tau=1,
        t0_range = None, # Time/10 if None
        alpha_min = 0.1,
        alpha_max = 1.0,
        alpha_step = 0.02,
        mulF0_min=1.0,
        mulF0_max=10.0
    ):
        self.Time = Time
        self.Time_shift = Time_shift
        self.tau = tau
        self.t0_range = t0_range or Time / 10.
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.alpha_step = alpha_step
        self.mulF0_min = mulF0_min
        self.mulF0_max = mulF0_max


class FRED(object):
    def __init__ (
        self,
        obs,
        sp,
        op=OptimizeParameters(),
        absolute=False,
        flux_obs='Lx',
        flux_model='Lx',
        flux_model_func=None,
        flux_fit_model='t0_Amplitude',
        cloptions={}
    ):
        self.obs = obs
        self.sp = sp
        self.op = op
        self.absolute = absolute
        self.flux_obs = flux_obs
        self.flux_model = flux_model
        self.flux_model_func = flux_model_func
        self.flux_fit_model = flux_fit_model
        self.cloptions = cloptions

        self.F0_0 = self.sp.Mdot_max * np.sqrt(G * Msun * self.sp.Mx * self.sp.r_out * Rsun) / (2. * np.pi)

        if 'err' in obs.dtype.names:
            # self.obs['err'] /= self.obs[self.flux_obs]
            pass
        else:
            self.obs = rec_append_fields( self.obs, 'err', np.ones_like(self.obs['DaP']) )
        if not self.absolute:
            self.obs[self.flux_obs] /= self.obs[self.flux_obs].max()
       
        plt.switch_backend('pdf')
        plt.yscale('log')

    def datadir( self, F0, alpha ):
        dirname = 'data/F_{F}/alpha_{alpha}'.format(F=F0, alpha=alpha)
        os.makedirs(dirname, exist_ok=True)
        return dirname
    
    def get_model( self, F0, alpha ):
        dirname = self.datadir(F0, alpha)
        model = np.genfromtxt( os.path.join(dirname, 'sum.dat'), names=True )
        if len( model.shape ) == 0:
            raise RuntimeError('Model file is empty')
        if model.shape[0] <= 3:
            raise RuntimeError('Need more dots in model file')
        if self.flux_model_func:
            model = rec_append_fields( model, self.flux_model+'_orig', model[self.flux_model] )
            model[self.flux_model] = self.flux_model_func( model[self.flux_model] )
        if not self.absolute:
            model[self.flux_model] /= model[self.flux_model].max()

        return model

    def fit_t0( self, F0, alpha, draw=False ):
        dirname = self.datadir(F0, alpha)

        subprocess.call(
            [
                './fred',
                '--alpha={}'.format(alpha),
                '--Mx={}'.format(self.sp.Mx),
                '--rout={}'.format(self.sp.r_out),
                #'--Mopt={}'.format(self.sp.Mopt),
                #'--period={}'.format(self.sp.Period),
                '--time={}'.format(self.op.Time),
                '--tau={}'.format(self.op.tau),
                '--dir={}'.format(dirname),
                '--F0={}'.format(F0),
                '--dilution={}'.format(self.sp.fcol),
            ] + [
                '--{}={}'.format(k,v) for k,v in self.cloptions.items() if not(v is None)
            ] + [
                '--{}'.format(k,v) for k, v in self.cloptions.items() if (v is None)
            ]
        )
        
        model = self.get_model(F0, alpha)
        if 'Mdot_orig' in model.dtype.names:
            Mdot0 = model['Mdot_orig'].max()
        else:
            Mdot0 = model['Mdot'].max()
        model_spline = interp1d( model['t'], model[self.flux_model], kind='cubic' )
        t0 = model['t'][ model[self.flux_model].argmax() ]
        i_min = self.obs['DaP'].searchsorted( -t0, side='left' )
        i_max = self.obs['DaP'].searchsorted( self.op.Time - t0 - self.op.t0_range, side='right' )
        obs_trunc = self.obs[i_min:i_max]

        Amplitude = 1.

        if self.flux_fit_model == 't0_Amplitude':
            try:
                params = curve_fit(
                    lambda t,t0,A: A * model_spline(t+t0),
                    obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'],
                    p0=(t0, Amplitude)
                )
                t0 = params[0][0]
                t0_err = params[1][0][0]
                Amplitude = params[0][1]
                Amplitude_err = params[1][1][0]
            except (ValueError, TypeError):
                t0_err = float('inf')
        elif self.flux_fit_model == 't0':
            try:
                params = curve_fit(
                    lambda t,t0: model_spline(t+t0),
                    obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'],
                    p0=(t0,)
                )
                t0 = params[0][0]
                t0_err = params[1][0][0]
            except (ValueError, TypeError):
                t0_err = float('inf')
        elif self.flux_fit_model == 'None' or self.flux_fit_model is None:
            obs_trunc = self.obs[ self.obs['DaP'] <= self.op.Time ]
            if model['t'][-1] < self.op.Time:
                t0_err = float('inf')
            else:
                t0_err = 0.
        else:
            raise RuntimeError('Wrong flux_fit_model')

        if draw:
            plt.cla()
            #plt.plot( obs_trunc['DaP'], obs_trunc[self.flux_obs], 'x' )
            plt.errorbar( self.obs['DaP'], self.obs[self.flux_obs], yerr=self.obs['err'], linestyle='None' )
            plt.plot( model['t']-t0, Amplitude * model[self.flux_model] )
            plt.savefig( os.path.join(dirname, 'curves.pdf'), format='pdf' )

        if t0_err == float('inf'):
            raise RuntimeError('Cannot fit with this parameters')
        
        return ( lambda t: Amplitude * model_spline(t+t0), model, obs_trunc, t0, Mdot0 )

    def chi2_fit_t0( self, F0, alpha ):
        print(F0/self.F0_0, alpha)
        try:
            model_spline, model, obs_trunc, t0, Mdot0 = self.fit_t0( F0, alpha )
        except RuntimeError:
            return float('inf')
        return chi2( model_spline, obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'] )

    def texp( self, model_spline, t0 ):
        return (self.op.Time-t0) / np.log( model_spline(0.) / model_spline(self.op.Time-t0) )

    def fit_alpha(self):
        #for F0 in 10**np.linspace(37.5, 39.5, 6):
        Mdot0 = self.sp.Mdot_max * 2
        F0 = self.sp.Mdot_max * np.sqrt(G * Msun * self.sp.Mx * self.sp.r_out * Rsun) / (2. * np.pi)
        alpha = 1. # Initial value for alpha
        while abs( Mdot0 / self.sp.Mdot_max - 1. ) > 1e-2:
            min_chi2 = float('inf')
            for alpha in np.arange( self.op.alpha_min, self.op.alpha_max + self.op.alpha_step, self.op.alpha_step ):
                alpha = np.around( alpha, decimals=int(np.ceil( np.log10(1/self.op.alpha_step) )) )
                try:
                    model_spline, model, obs_trunc, t0, Mdot0 = self.fit_t0( F0, alpha )
                except RuntimeError:
                    continue
                c = chi2( model_spline, obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'] )
                if c < min_chi2:
                    min_chi2 = c
                    alpha_min_chi2 = alpha
            alpha = alpha_min_chi2
            #minimize_result = minimize(
            #    lambda alpha: self.chi2_fit_t0(F0, alpha[0]),
            #    (alpha,),
            #    bounds=( (self.op.alpha_min, self.op.alpha_max) ),
            #    options={ 'maxiter': int( (self.op.alpha_max-self.op.alpha_min) // self.op.alpha_step ) },
            #)
            #if not minimize_result.success:
            #    raise RuntimeError("Error in minimize procedure")
            #alpha = minimize_result.x[0]
            #print(minimize_result)
            model_spline, model, obs_trunc, t0, Mdot0 = self.fit_t0( F0, alpha )
            H2R = interp1d( model['t'], model['H2R'], kind='cubic' )(t0)
            t_exp = self.texp( model_spline, t0 )
            print( F0, alpha, H2R, t0, t_exp, (H2R/0.05)**2 * alpha )
            F0 *= self.sp.Mdot_max / Mdot0

    def fit_F0(self, alpha):
        print(self.F0_0)
        minimize_result = differential_evolution(
            lambda x: self.chi2_fit_t0( x[0] * self.F0_0, alpha ),
            ( (self.op.mulF0_min, self.op.mulF0_max), ),
            tol=2e-2,
            popsize=7,
            maxiter=50,
            seed=np.random.mtrand.RandomState(),
        )
        print(minimize_result)
        F0 = self.F0_0 * minimize_result.x[0]
        return F0
    
    def fit_F0alpha(self):
        print(self.F0_0)
        minimize_result = differential_evolution(
            lambda x: self.chi2_fit_t0( x[0] * self.F0_0, x[1] ),
            (
                (self.op.mulF0_min, self.op.mulF0_max),
                (self.op.alpha_min, self.op.alpha_max)
            ),
            tol=2e-2,
            popsize=7,
            maxiter=200,
            seed=np.random.mtrand.RandomState(),
        )
        print( minimize_result )
        F0 = self.F0_0 * minimize_result.x[0]
        alpha = minimize_result.x[1]
        return( F0, alpha )
    
    def find_errors(self, F0, alpha, sigma_level=1):
        chi2_func = lambda x: self.chi2_fit_t0( x[0] * self.F0_0, x[1] )
        x0 = (F0/self.F0_0, alpha)
        bounds = (
            (self.op.mulF0_min, self.op.mulF0_max),
            (self.op.alpha_min, self.op.alpha_max)
        )
        dof = self.get_model(F0, alpha).shape[0]-2
        errors = chi2_errors( 
            chi2_func,
            x0, bounds,
            dof,
            sigma_level=sigma_level
        )
        # errors[0] *= self.F0_0
        return errors

    def print_params(self, F0, alpha, stream=sys.stdout, oneline=False, additional_fields={}):
        model_spline, model, obs_trunc, t0, Mdot0 = self.fit_t0(F0, alpha, draw=True)
        r_cold = self.sp.r_out
        r_hot_0 = r_cold * model[0]['Rhot2Rout']

        splines = {}
        for column in model.dtype.names[1:]:
            splines[column] = interp1d( model['t'], model[column], kind='cubic' )

        t0_eff = t0 if t0>0 else 0

        if oneline:
            stream.write(
                '{Mx:<6.1f} {Kerr:<5.3g} {Mopt:<4.1f} {Period:<7.4g} {r_out:<7.5f} {fullTime:<8g} {spectrum_fit:<25s} {t0range:<7g} {F0:<10.4e} {alpha:<7.5f} {initialcond:<11s} {ic_param:<8.1g} {Thot:<6g} {run_radius:<10d} {k_irr:<5.2g} {irr_enable:<10d} {fitdots:<7s} {fcol:<5.1f} {Chi2:<5.3g} {r_hot_0:<7.5f} {r_hot_max:<9.5g} {k_x_max:<9.4g} {H2r:<7.5f} {Mdot_max:<10.4e} {T0:<10.1f}\n'.format(
                    Mx = self.sp.Mx,
                    Kerr = self.cloptions.get('kerr') or 0.,
                    Mopt = self.sp.Mopt,
                    Period = self.sp.Period,
                    r_out = self.sp.r_out,
                    fullTime = self.op.Time,
                    spectrum_fit = additional_fields.get('spectrum_fit'),
                    t0range = self.op.t0_range,
                    F0 = F0,
                    alpha = alpha,
                    initialcond = self.cloptions.get('initialcond') or 'unknown',
                    ic_param = self.cloptions.get('powerorder') or float('inf'),
                    Thot = self.cloptions['Thot'],
                    run_radius = int( self.cloptions['Thot'] > 0. ),
                    k_irr = self.cloptions.get('Cirr'),
                    irr_enable = int( self.cloptions['Cirr'] > 0. ),
                    fitdots = self.flux_obs,
                    fcol = self.sp.fcol,
                    Chi2 = chi2( model_spline, obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'] ),
                    r_hot_0 = r_hot_0,
                    r_hot_max = float( splines['Rhot2Rout'](t0_eff) * r_cold ),
                    k_x_max = float( splines['kxout'](t0_eff) ),
                    H2r = float( splines['H2R'](t0_eff) ),
                    Mdot_max = Mdot0,
                    T0 = self.op.Time_shift - t0,
                )
            )
        else:
            stream.write('''
            F0 = {F0}
            alpha = {alpha}
            chi2 = {chi2}
            t0 = {t0}
            Mdot_max = {Mdot_max}
            r_cold = {r_cold}
            r_hot_0 = {r_hot_0}
            r_hot_max = {r_hot_max}
            k_x_max = {k_x_max}
            H2r = {H2r}
            dir = {path}
            '''.format(
                F0=F0,
                alpha=alpha,
                chi2=chi2( model_spline, obs_trunc['DaP'], obs_trunc[self.flux_obs], sigma=obs_trunc['err'] ),
                t0=t0,
                Mdot_max=Mdot0,
                r_cold=r_cold,
                r_hot_0 = r_hot_0,
                r_hot_max = splines['Rhot2Rout'](t0) * r_cold,
                k_x_max = splines['kxout'](t0),
                H2r = splines['H2R'](t0),
                path = self.datadir(F0, alpha),
            ))
