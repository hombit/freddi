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
* [Questions and comments](#questions-and-comments)
* [License](#license)

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. The
code is developed to simulate fast rise exponential decay (FRED) light curves of
low mass X-ray binaries (LMXBs) for the paper “Determination of the turbulent
parameter in the accretion disks: effects of self-irradiation in 4U 1543-47
during the 2002 outburst” by Lipunova & Malanchev (2016)
[arXiv:1610.01399](https://arxiv.org/abs/1610.01399).

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

-   [Boost](http://www.boost.org/) 1.56+

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
Freddi - numerical calculation of accretion disk evolution:

General options:
  -h [ --help ]                    Produce help message
  --prefix arg (=freddi)           Set prefix for output filenames. Output file
                                   with distribution of parameters over time is
                                   PREFIX.dat
  -d [ --dir ] arg (=.)            Choose the directory to write output files. 
                                   It should exist
  --fulldata                       Output files PREFIX_%d.dat with radial 
                                   structure for every time step. Default is to
                                   output only PREFIX.dat with global disk 
                                   parameters for every time step

Basic binary and disk parameters:
  -M [ --Mx ] arg (=5)             Mass of the central object, in the units of 
                                   solar masses
  --kerr arg (=0)                  Dimensionless Kerr parameter of the black 
                                   hole
  -a [ --alpha ] arg (=0.25)       Alpha parameter of Shakura-Sunyaev model
  --rin arg                        Inner radius of the disk, in the units of 
                                   the Schwarzschild radius of the central 
                                   object 2GM/c^2. If it isn't set then the 
                                   radius of ISCO orbit is used defined by --Mx
                                   and --kerr values
  --Mopt arg (=0.5)                Mass of the optical star, in units of solar 
                                   masses
  -P [ --period ] arg (=0.25)      Orbital period of the binary system, in 
                                   units of days
  -R [ --rout ] arg                Outer radius of the disk, in units of solar 
                                   radius. If it isn't set then the tidal 
                                   radius is used defined by --Mx, --Mopt and 
                                   --period values

Parameters of the disk model:
  -O [ --opacity ] arg (=Kramers)  Opacity law: Kramers (varkappa ~ rho / 
                                   T^7/2) or OPAL (varkappa ~ rho / T^5/2)
  --boundcond arg (=Teff)          Outer boundary movement condition
                                   
                                   Values:
                                     Teff: outer radius of the disk moves 
                                   inwards to keep photosphere temperature of 
                                   the disk larger than some value. This value 
                                   is specified by --Thot option
                                     Tirr: outer radius of the disk moves 
                                   inwards to keep irradiation flux of the disk
                                   larger than some value. The value of this 
                                   minimal irradiation flux is 
                                   [Stefan-Boltzmann constant] * Tirr^4, where 
                                   Tirr is specified by --Thot option
  --Thot arg (=0)                  Minimum photosphere or irradiation 
                                   temperature at the outer edge of the hot 
                                   disk, Kelvin. For details see --boundcond 
                                   description
  --F0 arg (=2e+38)                Initial maximum viscous torque in the disk, 
                                   dyn*cm. Can be overwritten via --Mdisk0 and 
                                   --Mdot0
  --Mdisk0 arg                     Initial disk mass, g. If both --F0 and 
                                   --Mdisk0 are specified then --Mdisk0 is 
                                   used. If both --Mdot0 and --Mdisk0 are 
                                   specified then --Mdot0 is used
  --Mdot0 arg                      Initial mass accretion rate through the 
                                   inner radius, g/s. If --F0, --Mdisk0 and 
                                   --Mdot0 are specified then --Mdot0 is used. 
                                   Works only when --initialcond is set to 
                                   sinusF or quasistat
  --initialcond arg (=powerF)      Type of the initial condition for viscous 
                                   torque F or surface density Sigma
                                   
                                   Values:
                                     powerF: F ~ xi^powerorder, powerorder is 
                                   specified by --powerorder option
                                     powerSigma: Sigma ~ xi^powerorder, 
                                   powerorder is specified by --powerorder 
                                   option
                                     sinusF: F ~ sin( xi * pi/2 )
                                     gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2),
                                   mu and sigma are specified by --gaussmu and 
                                   --gausssigma options
                                     quasistat: F ~ f(h/h_out) * xi * h_out/h, 
                                   where f is quasi-stationary solution found 
                                   in Lipunova & Shakura 2000. f(xi=0) = 0, 
                                   df/dxi(xi=1) = 0
                                   
                                   Here xi is (h - h_in) / (h_out - h_in)
                                   
  --powerorder arg (=6)            Parameter for the powerlaw initial condition
                                   distribution. This option works only with 
                                   --initialcond=powerF or powerSigma
  --gaussmu arg (=1)               Position of the maximum for Gauss 
                                   distribution, positive number not greater 
                                   than unity. This option works only with 
                                   --initialcond=gaussF
  --gausssigma arg (=0.25)         Width of for Gauss distribution. This option
                                   works only with --initialcond=gaussF

Parameters of self-irradiation:
  --Cirr arg (=0)                  Irradiation factor
  --irrfactortype arg (=const)     Type of irradiation factor Cirr
                                   
                                   Values:
                                     const: doesn't depend on disk shape:
                                   [rad. flux] = Cirr  L / (4 pi r^2)
                                     square: Cirr depends on the disk relative 
                                   half-thickness:
                                   [rad. flux] = Cirr (z/r)^2 L / (4 pi r^2)
                                   
                                   Here L is bolometric Luminosity:
                                   L = eta Mdot c^2

Parameters of optical magnitudes calculation:
  --colourfactor arg (=1.7)        Colour factor to calculate X-ray flux
  --emin arg (=1)                  Minimum energy of X-ray band, keV
  --emax arg (=12)                 Maximum energy of X-ray band, keV
  -i [ --inclination ] arg (=0)    Inclination of the system, degrees
  --distance arg (=10)             Distance to the system, kpc
  --lambda arg                     Wavelength to calculate Fnu, Angstrom. You 
                                   can use this option multiple times. For each
                                   lambda one additional column with values of 
                                   spectral flux density Fnu [erg/s/cm^2/Hz] is
                                   produced

Parameters of disk evolution calculation:
  -T [ --time ] arg (=50)          Time interval to calculate evolution, days
  --tau arg (=0.25)                Time step, days
  --Nx arg (=1000)                 Size of calculation grid
  --gridscale arg (=log)           Type of grid for angular momentum h: log or 
                                   linear

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
  --kerr=0.4 --Cirr=0.001 --opacity=OPAL --initialcond=quasistat
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
of the black using the approximation by [Suleimanov et al
(2008)](http://adsabs.harvard.edu/abs/2008A%26A...491..267S).

The parameters at the disk central plane are defined by the analytic
approximations ([Suleimanov et al.
2007](http://adsabs.harvard.edu/abs/2007ARep...51..549S)), valid for the
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

Questions and comments
----------------------

If you have any problems, questions, or comments, please address them to
[Issues](https://github.com/hombit/freddi/issues) or to hombit\@gmail.com

License
-------

Copyright (c) 2016, Konstantin L. Malanchev & Galina V. Lipunova.

`Freddi` is distributed under the terms of the
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2016) [arXiv:1610.01399](https://arxiv.org/abs/1610.01399)
