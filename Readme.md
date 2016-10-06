`Freddi` — compute FRED-like light curves of LMXB
=================================================

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
further installation instructions and go straight to the Usage section, using
following string instead of `./freddi`.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
docker run -v "`pwd`":/data --rm -ti hombit/freddi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Requirements

-   [Boost.Program\_options](http://www.boost.org/doc/libs/release/doc/html/program_options.html)

-   C++ compiler with C++11 support, e.g. `gcc` version 4.8+ or `clang` 3.4+

`Freddi` was tested on Linux and macOS but it should work on Windows as well.

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

`Freddi` runs from the command line with optionally set arguments decribed
below.`Freddi` always outputs file `freddi.dat` with  distribution of various
physical values over time. If `--fulldata` is specified then files
`freddi_%d.dat` for each time step are created in the same directory with
snapshot radial distributions. These `*.dat` data-files contain
whitespace-separated data columns with header lines starting with `#` symbol.
You can use another prefix instead of `freddi` with `--prefix` option and change
output directory with `--dir` option. If you choose Docker way and would like to
specify the directory, then avoid to use `--dir` option and just replace ``
"`pwd`" `` with some local path (for more details see [Docker
documentation](https://docs.docker.com/engine/tutorials/dockervolumes)).

See the full list of command line options with `--help` option. Default values
are given in brackets.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ ./freddi --help
Freddi - numerical calculation of accretion disk evolution:

General options:
  -h [ --help ]                         Produce help message
  --prefix arg (=freddi)                Prefix for output filenames. Output 
                                        file with distribution of parameters 
                                        over time is PREFIX.dat
  -d [ --dir ] arg (=.)                 Directory to write output files. It 
                                        should exist
  --fulldata                            Output files PREFIX_%d.dat with radial 
                                        structure for every time step. Default 
                                        is to output only PREFIX.dat with 
                                        global disk parameters for every time 
                                        step

Basic binary and disk parameters:
  -M [ --Mx ] arg (=10)                 Mass of the central object, in the 
                                        units of solar masses
  --kerr arg (=0)                       Dimensionless Kerr parameter of the 
                                        black hole
  -a [ --alpha ] arg (=0.25)            Alpha parameter of Shakura-Sunyaev 
                                        model
  --rin arg                             Inner radius of the disk, in the units 
                                        of the Schwarzschild radius of the 
                                        central object 2GM/c^2. If it isn't set
                                        then the radius of ISCO orbit is used 
                                        defined by --Mx and --kerr values
  --Mopt arg (=1)                       Mass of the optical star, in units of 
                                        solar masses
  -P [ --period ] arg (=1)              Orbital period of the binary system, in
                                        units of days
  -R [ --rout ] arg (=4.3302577068820627)
                                        Outer radius of the disk, in units of 
                                        solar radius. If it isn't set then the 
                                        tidal radius is used defined by --Mx, 
                                        --Mopt and --period values
  -i [ --inclination ] arg (=0)         Inclination of the system, degrees

Parameters of the disk model:
  -O [ --opacity ] arg (=Kramers)       Opacity law: Kramers (varkappa ~ rho / 
                                        T^7/2) or OPAL (varkappa ~ rho / T^5/2)
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
  --F0 arg (=1e+36)                     Initial viscous torque at the outer 
                                        boundary of the disk, dyn*cm
  --Mdot0 arg (=0)                      Initial mass accretion rate through the
                                        inner radius, g/s. If both --F0 and 
                                        --Mdot0 are specified then --Mdot0 is 
                                        used. Works only when --initialcond is 
                                        set to sinusF or quasistat
  --initialcond arg (=power)            Type of the initial condition for 
                                        viscous torque F or surface density 
                                        Sigma
                                        
                                        Values:
                                          powerF: F ~ xi^powerorder, powerorder
                                        is specified by --powerorder option
                                          powerSigma: Sigma ~ xi^powerorder, 
                                        powerorder is specified by --powerorder
                                        option
                                          sinusF: F ~ sin( xi * pi/2 )
                                          quasistat: F ~ f(h/h_out) * xi * 
                                        h_out/h, where f is quasi-stationary 
                                        solution found in Lipunova & Shakura 
                                        2000. f(xi=0) = 0, df/dxi(xi=1) = 0
                                        
                                        Here xi is (h - h_in) / (h_out - h_in)
                                        
  --powerorder arg (=6)                 Parameter for the powerlaw initial 
                                        condition distribution. This option 
                                        works only with --initialcond=powerF or
                                        powerSigma

Parameters of X-ray emission:
  --Cirr arg (=0)                       Irradiation factor
  --irrfactortype arg (=const)          Type of irradiation factor Cirr
                                        
                                        Values:
                                          const: doesn't depend on disk shape:
                                        [rad. flux] = Cirr  L / (4 pi r^2)
                                          square: disk has polynomial shape:
                                        [rad. flux] = Cirr (z/r)^2 L / (4 pi 
                                        r^2)
                                        
  --colourfactor arg (=1.7)             Colour factor
  --emin arg (=1)                       Lower bound of X-ray band, keV
  --emax arg (=12)                      Upper bound of X-ray band, keV

Parameters of optical magnitudes calculation:
  --distance arg (=10)                  Distance to the system, kpc

Parameters of disk evolution calculation:
  -T [ --time ] arg (=25)               Time interval to calculate evolution, 
                                        days
  --tau arg (=0.25)                     Time step, days
  --Nx arg (=1000)                      Size of calculation grid
  --gridscale arg (=log)                Type of grid for angular momentum h: 
                                        log or linear
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
be found on solving the equation. The accretion rate defined the X-ray
luminosity of the source.

Standard model for the accretion disk is implied, which is developed by [Shakura
& Sunyaev (1973)](http://adsabs.harvard.edu/abs/1973A%26A....24..337S). The
metric is Newtonian which is accurate enough for the problem. The boundary
conditions in the disk are the zero stress at the inner boundary and the zero
accretion rate at the outer boundary. The conditions are suitable for the
outbursts in X-ray binary transients with black holes. The vertical structure of
the disk is solved in the code and valid for the effective surface temperatures
from 10 000 to 100 000 K, approximately. It is assumed that the gas pressure
dominates in the disk, the gas is completely ionized, and the photon opacity is
defined by free-free and free-bound transitions (Kramers' type of
OPAL-approximated opacity for chemical composition with solar abundances).

The outer radius of the disk can be stationary or moving in the approximation
that the evolution goes through the quasi-stationary states. Different options
are implemented to control the position of the outer radius.

The initial distribution of the matter in the disk should be specified with
`--initialcond` option. `Freddi` can start from several types of initial
distributions: power-law distribution of the surface density
`--initialcond=powerSigma`, power-law `--initialcond=powerF` or sinus-law
`--initialcond=sinusF` distribution of the viscous torque, quasi-stationary
distribution `--initialcond=quasistat`. The choice of the initial distribution
defines what type of evolution is to be calculated. Starting from the
quasi-stationary distribution, the solution describes the decaying part of the
outburst, otherwise, the rise to the peak is also computed.

Fitting parameters of the disks in X-ray transients is one of the `Freddi`
goals.

License
-------

Copyright (c) 2016, Konstantin L. Malanchev & Galina V. Lipunova.

`Freddi` is distributed under the terms of the
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2016) [arXiv:1610.01399](https://arxiv.org/abs/1610.01399)
