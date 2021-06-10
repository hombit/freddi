`Freddi` — compute FRED-like light curves of LMXB
=================================================

Table of contents
-----------------

* [Overview](#overview)
* [Installation](#installation)
  * [Docker](#docker)
  * [Requirements](#requirements)
  * [Get and compile source files](#get-and-compile-source-files)
* [Usage](#usage)
  * [Options](#options)
  * [Output values](#output-values)
  * [Example](#example)
* [Physical Background](#physical-background) 
* [Accretion disk wind](#accretion-disk-wind) 
  * [Compton-heated wind](#compton-heated-wind)
* [Questions and comments](#questions-and-comments)
* [License](#license)
* [BibTex](#bibtex)

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. The
code is developed to simulate fast-rise exponential-decay (FRED) light curves of
low mass X-ray binaries (LMXBs) for the paper “Determination of the turbulent
parameter in the accretion disks: effects of self-irradiation in 4U 1543-47
during the 2002 outburst” by Lipunova & Malanchev (2017)
[2017MNRAS.468.4735L](http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L).

Installation
------------

### Docker

If you are familiar with [Docker](http://docker.com) then you can skip all
further installation instructions and go straight to the [Usage section](#usage), using
following string instead of `./freddi`.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run -v "`pwd`":/data --rm -ti hombit/freddi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Requirements

-   [Boost](http://www.boost.org/) 1.57+

-   [Make](https://en.wikipedia.org/wiki/Make_(software))

-   C++ compiler with C++11 support, e.g. `gcc` version 4.8+ or `clang` 3.4+

Requirements installation on Debian based systems (e.g. Ubuntu):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
apt-get install g++ make libboost-all-dev
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Red Hat based systems (e.g. Fedora):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dnf install gcc-c++ make boost-devel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Freddi` was tested on Linux and macOS.

### Get and compile source files

First go to the path where `Freddi` directory will be located. Then download and
compile it:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
git clone https://github.com/hombit/freddi.git
cd freddi
make
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now it should be executable file `./freddi` in the current directory. If you’d
like to install it to `/usr/local/bin` then do

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sudo make install
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage
-----

`Freddi` runs from the command line with inline options and/or with a configuration file. `Freddi`
outputs file `freddi.dat` with distribution of various physical values over
time. If `--fulldata` is specified then files `freddi_%d.dat` for each time step
are created in the same directory with snapshot radial distributions. These
data-files contain whitespace-separated data columns with header lines starting
with `#` symbol. You can set another prefix instead of `freddi` with `--prefix`
option and change the output directory with `--dir` option. If you choose the
Docker way and would like to specify the directory, then avoid using `--dir`
option and just replace `` "`pwd`" `` with some local path (for more details see
[Docker documentation](https://docs.docker.com/engine/tutorials/dockervolumes)).

### Options

The full list of command line options is viewed with `--help` option. Default
values are given in brackets.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ ./freddi --help
Freddi: numerical calculation of accretion disk evolution:

General options:
  -h [ --help ]                         Produce help message
  --config arg                          Set additional configuration filepath
  --prefix arg (=freddi)                Set prefix for output filenames. Output
                                        file with distribution of parameters 
                                        over time is PREFIX.dat
  -d [ --dir ] arg (=.)                 Choose the directory to write output 
                                        files. It should exist
  --precision arg (=12)                 Number of digits to print into output 
                                        files
  --fulldata                            Output files PREFIX_%d.dat with radial 
                                        structure for every time step. Default 
                                        is to output only PREFIX.dat with 
                                        global disk parameters for every time 
                                        step

Basic binary and disk parameter:
  -a [ --alpha ] arg                    Alpha parameter of Shakura-Sunyaev 
                                        model
  --alphacold arg                       Alpha parameter of cold disk, currently
                                        it is used only for Sigma_minus, see 
                                        --Qirr2Qvishot. Default is --alpha 
                                        values divided by ten
  -M [ --Mx ] arg                       Mass of the central object, in the 
                                        units of solar masses
  --kerr arg (=0)                       Dimensionless Kerr parameter of the 
                                        black hole
  --Mopt arg                            Mass of the optical star, in units of 
                                        solar masses
  --rochelobefill arg (=1)              Dimensionless factor describing a size 
                                        of the optical star. Polar radius of 
                                        the star is rochelobefill * (polar 
                                        radius of critical Roche lobe)
  --Topt arg (=0)                       Thermal temperature of the optical 
                                        star, in units of kelvins
  -P [ --period ] arg                   Orbital period of the binary system, in
                                        units of days
  --rin arg                             Inner radius of the disk, in the units 
                                        of the gravitational radius of the 
                                        central object GM/c^2. If it isn't set 
                                        then the radius of ISCO orbit is used 
                                        defined by --Mx and --kerr values
  -R [ --rout ] arg                     Outer radius of the disk, in units of 
                                        solar radius. If it isn't set then the 
                                        tidal radius is used defined by --Mx, 
                                        --Mopt and --period values
  --risco arg                           Innermost stable circular orbit, in 
                                        units of gravitational radius of the 
                                        central object GM/c^2. If it isn't set 
                                        then the radius of ISCO orbit is used 
                                        defined by --Mx and --kerr values

Parameters of the disk mode:
  -O [ --opacity ] arg (=Kramers)       Opacity law: Kramers (varkappa ~ rho / 
                                        T^7/2) or OPAL (varkappa ~ rho / T^5/2)
  --Mdotout arg (=0)                    Accretion rate onto the disk through 
                                        its outer radius
  --boundcond arg (=Teff)               Outer boundary movement condition
                                        
                                        Values:
                                          Teff: outer radius of the disk moves 
                                        inwards to keep photosphere temperature
                                        of the disk larger than some value. 
                                        This value is specified by --Thot 
                                        option
                                          Tirr: outer radius of the disk moves 
                                        inwards to keep irradiation flux of the
                                        disk larger than some value. The value 
                                        of this minimal irradiation flux is 
                                        [Stefan-Boltzmann constant] * Tirr^4, 
                                        where Tirr is specified by --Thot 
                                        option
  --Thot arg (=0)                       Minimum photosphere or irradiation 
                                        temperature at the outer edge of the 
                                        hot disk, Kelvin. For details see 
                                        --boundcond description
  --Qirr2Qvishot arg (=0)               Minimum Qirr / Qvis ratio at the outer 
                                        edge of the hot disk to switch 
                                        evolution from temperature-based regime
                                        to Sigma_minus-based regime (see Eq. 
                                        A.1 in Lasota et al. 2008, --alphacold 
                                        value is used as alpha parameter)
  --initialcond arg (=powerF)           Type of the initial condition for 
                                        viscous torque F or surface density 
                                        Sigma
                                        
                                        Values:
                                          powerF: F ~ xi^powerorder, powerorder
                                        is specified by --powerorder option
                                          linearF: F ~ xi, specific case of 
                                        powerF but can be normalised by 
                                        --Mdot0, see its description for 
                                        details  powerSigma: Sigma ~ 
                                        xi^powerorder, powerorder is specified 
                                        by --powerorder option
                                          sineF: F ~ sin( xi * pi/2 )
                                          gaussF: F ~ exp(-(xi-mu)**2 / 2 
                                        sigma**2), mu and sigma are specified 
                                        by --gaussmu and --gausssigma options
                                          quasistat: F ~ f(h/h_out) * xi * 
                                        h_out/h, where f is quasi-stationary 
                                        solution found in Lipunova & Shakura 
                                        2000. f(xi=0) = 0, df/dxi(xi=1) = 0
                                        
                                        Here xi is (h - h_in) / (h_out - h_in)
                                        
  --F0 arg                              Initial maximum viscous torque in the 
                                        disk, dyn*cm. Can be overwritten via 
                                        --Mdisk0 and --Mdot0
  --Mdisk0 arg                          Initial disk mass, g. If both --F0 and 
                                        --Mdisk0 are specified then --Mdisk0 is
                                        used. If both --Mdot0 and --Mdisk0 are 
                                        specified then --Mdot0 is used
  --Mdot0 arg                           Initial mass accretion rate through the
                                        inner radius, g/s. If --F0, --Mdisk0 
                                        and --Mdot0 are specified then --Mdot0 
                                        is used. Works only when --initialcond 
                                        is set to linearF, sinusF or quasistat
  --powerorder arg                      Parameter for the powerlaw initial 
                                        condition distribution. This option 
                                        works only with --initialcond=powerF or
                                        powerSigma
  --gaussmu arg                         Position of the maximum for Gauss 
                                        distribution, positive number not 
                                        greater than unity. This option works 
                                        only with --initialcond=gaussF
  --gausssigma arg                      Width of for Gauss distribution. This 
                                        option works only with 
                                        --initialcond=gaussF

Parameters of self-irradiation.
Qirr = Cirr * (H/r / 0.05)^irrindex * L * psi / (4 pi R^2), where psi is angular distrbution of X-ray radiation:
  --Cirr arg (=0)                       Irradiation factor for the hot disk
  --irrindex arg (=0)                   Irradiation index for the hot disk
  --Cirrcold arg (=0)                   Irradiation factor for the cold disk
  --irrindexcold arg (=0)               Irradiation index for the cold disk
  --h2rcold arg (=0.050000000000000003) Seme-height to radius ratio for the 
                                        cold disk
  --angulardistdisk arg (=plane)        Angular distribution of the disk X-ray 
                                        radiation. Values: isotropic, plane

Parameters of flux calculation:
  --colourfactor arg (=1.7)             Colour factor to calculate X-ray flux
  --emin arg (=1)                       Minimum energy of X-ray band, keV
  --emax arg (=12)                      Maximum energy of X-ray band, keV
  --staralbedo arg (=0)                 Part of X-ray radiation reflected by 
                                        optical star, (1 - albedo) heats star's
                                        photosphere. Used only when --starflux 
                                        is specified
  -i [ --inclination ] arg (=0)         Inclination of the system, degrees
  --ephemerist0 arg (=0)                Ephemeris for the time of the minimum 
                                        of the orbital light curve T0, phase 
                                        zero corresponds to inferior 
                                        conjunction of the optical star, days
  --distance arg                        Distance to the system, kpc
  --colddiskflux                        Add Fnu for cold disk into output file.
                                        Default output is for hot disk only
  --starflux                            Add Fnu for optical star into output 
                                        file. Mx, Mopt and period must be 
                                        specified, see also Topt and starlod 
                                        options. Default output is for hot disk
                                        only
  --lambda arg                          Wavelength to calculate Fnu, Angstrom. 
                                        You can use this option multiple times.
                                        For each lambda one additional column 
                                        with values of spectral flux density 
                                        Fnu [erg/s/cm^2/Hz] is produced
  --passband arg                        Path of a file containing tabulated 
                                        passband, the first column for 
                                        wavelength in Angstrom, the second 
                                        column for transmission factor, columns
                                        should be separated by spaces

Parameters of disk evolution calculation:
  --inittime arg (=0)                   Initial time moment, days
  -T [ --time ] arg                     Time interval to calculate evolution, 
                                        days
  --tau arg                             Time step, days
  --Nx arg (=1000)                      Size of calculation grid
  --gridscale arg (=log)                Type of grid for angular momentum h: 
                                        log or linear
  --starlod arg (=3)                    Level of detail of the optical star 3-D
                                        model. The optical star is represented 
                                        by a triangular tile, the number of 
                                        tiles is 20 * 4^starlod

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Also you can use `freddi.ini` configuration file to store options. This [INI
file](https://en.wikipedia.org/wiki/INI_file) contains lines `option=value`,
option names are the as provided by the help message above. Command line option
overwrites configuration file option. For example, [see
default](https://github.com/hombit/freddi/blob/master/freddi.ini) `freddi.ini`.

Paths where this file is searched are `./freddi.ini` (execution path),
`$HOME/.config/freddi/freddi.ini`, `/usr/local/etc/freddi.ini` and
`/etc/freddi.ini`. You can provide configuration file to Docker container as a
volume: `` -v "`pwd`/freddi.ini":/etc/freddi.ini ``.

### Output values

`Freddi` outputs time; the accretion rate; the mass of the hot part of the disk;
the outer radius of the hot zone; the irradiation factor; the relative
half-height, effective and irradiation temperature, ratio of the irradiation to
viscous flux at the outer radius of the hot zone; X-ray luminosity (erg/s) in
the band from E\_min to E\_max (`--emin` and `--emax` options); the optical
magnitudes in *U*, *B*, *V*, *R*, *I*, and *J* band ([Allen's Astrophysical
Quantities, Cox 2015](http://www.springer.com/book/9780387951898)); the spectral density flux (erg/s/cm^2/Hz) at some wavelengths set by one or more `--lambda` options.

Snapshot distributions at each time step, if produced, contain the following
data: radial coordinate in terms of the specific angular momentum, radius,
viscous torque, surface density, effective temperature Teff, viscous temperature
Tvis, irradiation temperature Tirr, and the absolute half-height of the disk.

### Example

The following arguments instruct `Freddi` to calculate the decay of the outburst
in the disk with the constant outer radius equal to 1 solar radius. The Kerr
black hole at the distance of 5 kpc has the mass of 9 solar masses, and the
Kerr's parameter is 0.4. The outer disk is irradiated with Cirr=1e-3.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
./freddi --alpha=0.5 --Mx=9 --rout=1 --time=50 --tau=0.25 --dir=data/ \
  --F0=2e+37 --colourfactor=1.7 --Nx=1000 --distance=5 --gridscale=log \
  --kerr=0.4 --Cirr=0.001 --opacity=OPAL --initialcond=quasistat \
  --wind=Woods1996 --Xi_max=10 --Tic=1e8 --M_pow=1 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Physical Background
-------------------

`Freddi` — Fast Rise Exponential Decay: accretion Disk model Implementation — is
designed to solve the differential equation of the viscous evolution of the
Shakura-Sunyaev accretion disk in a stellar binary system. Shakura-Sunyaev disk
is the standard model of accretion of plasma onto the cosmic bodies, like
neutron stars or black holes. Viscous evolution of the accretion disks exibits
itself, for example, in X-ray outbursts of binary stars. Usually, the outbursts
last for several tens of days and many of them are observed by orbital
observatories.

The basic equation of the viscous evolution relates the surface density and
viscous stresses and is of diffusion type. Evolution of the accretion rate can
be found on solving the equation. The distribution of viscous stresss defines
the emission from the source.

The standard model for the accretion disk is implied, which is developed by
[Shakura & Sunyaev (1973)](http://adsabs.harvard.edu/abs/1973A%26A....24..337S).
The inner boundary of the disk is at the ISCO or can be explicitely set. The
boundary conditions in the disk are the zero stress at the inner boundary and
the zero accretion rate at the outer boundary. The conditions are suitable
during the outbursts in X-ray binary transients with black holes.

In a binary system, the accretion disk is radially confined. In `Freddi`, the
outer radius of the disk can be set explicitely or calculated as the position of
the tidal truncation radius following [Paczynski
(1997)](http://adsabs.harvard.edu/abs/1977ApJ...216..822P) for small mass ratios
of the black using the approximation by [Suleimanov et al. (2008)](http://adsabs.harvard.edu/abs/2008A%26A...491..267S).

The parameters at the disk central plane are defined by the analytic
approximations ([Suleimanov et al. 2007](http://adsabs.harvard.edu/abs/2007ARep...51..549S)), valid for the
effective surface temperatures from 10 000 to 100 000 K, approximately. It is
assumed that the gas pressure dominates, the gas is completely ionized, and the
photon opacity is defined by the free-free and free-bound transitions. Opacity
law is for the solar element abundancies and can be chosen from two types: (1)
Kramers' opacity: kappa = 5e24 rho/T\^(7/2) cm2/g (2) approximation to OPAL
tables: kappa = 1.5e20 rho/T\^(5/2) cm2/g ([Bell & Lin
1994](http://adsabs.harvard.edu/abs/1994ApJ...427..987B))

The disk at each radius is in the "hot" state if the gas is completely ionized.
Otherwise, the disk is considered to be "cold" locally. Alpha-parameter in the
cold parts of the disk is appreciably lower than in the hot parts. Thus the
viscous evolution of the disk should proceed more effectively in the hot parts
of the disk. To simulate this, `Freddi` has an option to control the outer
radius of the hot evolving disk. We assume that the evolution goes through the
quasi-stationary states in the hot zone of variable size. By default, the hot
zone has the constant size, equal to the tidal radius.

The initial distribution of the matter in the disk should be specified with
`--initialcond` option. `Freddi` can start from several types of initial
distributions: power-law distribution of the surface density
`--initialcond=powerSigma`, power-law `--initialcond=powerF` or sinus-law
`--initialcond=sinusF` distribution of the viscous torque, quasi-stationary
distribution `--initialcond=quasistat`. The choice of the initial distribution
defines what type of evolution is to be calculated.

Starting from the quasi-stationary or `sinusF` distribution, the solution
describes the decaying part of the outburst. Zero-time accretion rate through
the inner edge can be set. In other cases, the rise to the peak is also
computed. Then, initial value of viscous torque at the outer radius (can be set
by `--F0`) defines uniquely the initial mass of the disk.

Self-irradiation by the central X-rays heats the outer parts of the disk. A
fraction of the bolometric flux is supposed to illuminate the disk surface. This
results in the larger size of the hot disk if such model is assumed. Also, the
optical flux is increased because the flux outgoing from the disk surface is
proportional to Teff\^4 = Tvis\^4+Tirr\^4. Self-irradiation of the disk is
included in the computation if irradiation parameter is not zero. The simplest
way is to set a constant irradiation factor `--Cirr` (the studies of X-ray novae
suggest the range for Cirr 1e-5—5e-3).

Observed flux depends on the distance to the source and the inclination of the
disk plane. The inclination angle is the angle between the line of sight and the
normal to the disk. The flux, emitted from the disk surface, is defined by the
sum of the visous and irradiating flux, where the viscous flux is calculated
taking into account general relativity effects near the black hole, following
[Page & Thorne
(1974)](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1974ApJ...191..499P)
and [Riffert & Herold (1995)](http://adsabs.harvard.edu/cgi-bin/nph-bib\_query?bibcode=1995ApJ...450..508R).


Accretion disk wind
-------------------

Presumably, during an outburst there is an outflow in the form of a wind from the 
accretion disk around the compact object. The presence of such a wind in the LMXBs is supported
by modern observations indicating the expansion of ionized matter. Such an outflow of matter,
being an additional  source of angular momentum transfer in the disk, can strongly influence 
its viscous evolution. 

However, the nature of such winds and their  physical characteristics are an open question.
Namely, there are three mechanisms which are considered:
heating of matter by central radiation in optically thin regions of the disk 
([Begelman et al. 1983](https://ui.adsabs.harvard.edu/abs/1983ApJ...271...70B), 
[Shields et al. 1986](https://ui.adsabs.harvard.edu/abs/1986ApJ...306...90S), 
[Woods et al. 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...461..767W)), 
the pressure of the magnetic field of the disk
([Blandford & Payne 1982](https://ui.adsabs.harvard.edu/abs/1982MNRAS.199..883B), 
[Habibi & Abbassi 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...887..256H), 
[Nixon & Pringle 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...628A.121N))
and the pressure of local radiation at super-Eddington accretion rates
([Shakura & Sunyaev 1973](http://adsabs.harvard.edu/abs/1973A%26A....24..337S), 
[Proga & Kallman 2002](https://ui.adsabs.harvard.edu/abs/2002ApJ...565..455P)). 

`Freddi` is modernized in such a way that it is able to solve the viscous evolution 
equation with an inhomogeneous term that is responsible for the presence of the disk wind.
This term is the dependence of the surface density of the wind mass-loss rate on 
the distance along the disk's surface. Different forms of such dependence correspond 
to different wind models, and to different classes within `Freddi`. 


One can choose a wind model by setting the
`--wind` option. The thermal wind model ([Woods et al. 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...461..767W)),
which implies that the outflow of matter occurs due to the heating of the outer parts of the disk
by a central radiation source, can be chosen by setting `--wind=Woods1996`. The option `--wind=Janiuk15`
corresponds to the model from work [Janiuk et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..92J)
where the wind is started in the super-Eddington regime. You can also select the `--wind=ToyWind` option, 
which corresponds to a toy wind model when the user sets the wind strength relatively to the accretion rate using the option `--windpow`.

### Compton-heated wind

At the moment, `Freddi` is more focused on simulating outbursts taking into account the thermal wind (`--wind=Woods1996` option). 
For a better understanding, let's discuss a little the physics of the process of launching such a wind 
and its parameters in the code.

In the standard accretion disk model by [Shakura & Sunyaev (1973)](http://adsabs.harvard.edu/abs/1973A%26A....24..337S) 
the disk is concave, and, as a result, the disk surface is exposed to the central radiation, 
which heats the disk material. As a result, the heated matter, starting from a certain radius, 
begins to leave the accretion disk. This process of heating the matter of the accretion disk by means of Compoton
processes was developed in [Begelman et al. (1983)](https://ui.adsabs.harvard.edu/abs/1983ApJ...271...70B) and 
[Shields et al. (1986)](https://ui.adsabs.harvard.edu/abs/1986ApJ...306...90S). 
In a later work [Woods et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...461..767W), 
two-dimensional magnetohydrodynamic calculations were performed and the 
results of [Shields et al. (1986)](https://ui.adsabs.harvard.edu/abs/1986ApJ...306...90S) were generalized. 
[Woods et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...461..767W) give an expression for the surface density of the mass 
loss rate as a function of distance along the disk's surface. This function is used in `Freddi` 
to taking thermal wind into account.

Choosing option `--wind=Woods1996`, it is necessary to set the value of the ionization parameter Xi
(which is proportional to the ratio of the radiation and gas pressures) by the option `--Xi_max` and the Compoton temperature T_IC 
(which determines the hardness of the irradiating spectrum and the size of the region where the wind operates) by the option `--Tic`. 


Questions and comments
----------------------

If you have any problems, questions, or comments, please address them to
[Issues](https://github.com/hombit/freddi/issues) or to hombit\@gmail.com

License
-------

Copyright (c) 2016–2021, Konstantin L. Malanchev, Galina V. Lipunova & Artur L. Avakyan.

`Freddi` is distributed under the terms of the
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2017) [2017MNRAS.468.4735L](http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L)

BibTex
------
```bibtex
@ARTICLE{2017MNRAS.468.4735L,
   author = {{Lipunova}, G.~V. and {Malanchev}, K.~L.},
    title = "{Determination of the turbulent parameter in accretion discs: effects of self-irradiation in 4U 1543{\minus}47 during the 2002 outburst}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1610.01399},
 primaryClass = "astro-ph.HE",
 keywords = {accretion, accretion discs, methods: numerical, binaries: close, stars: black holes, X-rays: individual: 4U 1543-47},
     year = 2017,
    month = jul,
   volume = 468,
    pages = {4735-4747},
      doi = {10.1093/mnras/stx768},
   adsurl = {http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
