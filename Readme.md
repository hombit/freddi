`Freddi` — compute FRED-like light curves of LMXB
=================================================

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. This
code is developed to simulate fast rise exponential decay (FRED) light curves of
low mass X-ray binaries (LMXBs) for the paper Lipunova & Malanchev (2016) (in
prep.).

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

`Freddi` always outputs `freddi.dat` file with temporal distribution of various
physical values. If `--fulldata` is specified then `freddi_%d.dat` files for
each time step are outputted to the same directory with various radial
distributions. These `*.dat` data-files contain whitespace-seporated data
columns with header lines started with `#` symbol. You can use another prefix
instead of `freddi` with `--prefix` option and change output directory with
`--dir` option. If you chose Docker way and would like to specify directory then
avoid to use `--dir` option and just replace `` "`pwd`" `` with some local path
(for more details see [Docker
documentation](https://docs.docker.com/engine/tutorials/dockervolumes)).

See full list of command line options with `--help` option.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ ./freddi --help
Freddi - numerical calculation of accretion disc evolution:

General options:
  -h [ --help ]                         Produce help message
  --prefix arg (=freddi)                Prefix for output filenames. File with 
                                        temporal distributions of parameters is
                                        PREFIX.dat
  -d [ --dir ] arg (=.)                 Directory to write output files. It 
                                        should exist
  --fulldata                            Output files PREFIX_%d.dat with radial 
                                        structure for every computed time step.
                                        Default is to output only PREFIX.dat 
                                        with global disk parameters for every 
                                        time step

Basic binary and disc parameters:
  -M [ --Mx ] arg (=10)                 Mass of the central object, solar 
                                        masses
  --kerr arg (=0)                       Kerr parameter of the black hole
  -a [ --alpha ] arg (=0.25)            Alpha parameter
  --rin arg                             Internal radius of the disk, 
                                        Schwarzschild radii of the central 
                                        object. If it isn't setted then it will
                                        be calculated as radius of ISCO orbit 
                                        using --Mx and --kerr values
  --Mopt arg (=1)                       Mass of optical star, solar masses
  -P [ --period ] arg (=1)              Orbital period of binary system, days
  -R [ --rout ] arg (=4.3302577068820627)
                                        Outer radius of the disk, solar radii. 
                                        If it isn't setted then it will be 
                                        calculated as tidal radius using --Mx, 
                                        --Mopt and --period
  -i [ --inclination ] arg (=0)         Inclination of the system, degrees

Parameters of the disc model:
  -O [ --opacity ] arg (=Kramers)       Opacity law: Kramers (varkappa ~ rho / 
                                        T^7/2) or OPAL (varkappa ~ rho / T^5/2)
  --boundcond arg (=Teff)               Outer boundary movement condition
                                        
                                        Values:
                                          Teff: outer radius of the disc moves 
                                        inside to keep photosphere temperature 
                                        of the disc larger than some value. 
                                        This value is specified by --Thot 
                                        option
                                          Tirr: outer radius of the disc moves 
                                        inside to keep irradiation flux of the 
                                        disc larger than some value. The value 
                                        of this minimal irradiation flux is 
                                        [Stefan-Boltzmann constant] * Tirr^4, 
                                        where Tirr is specified by --Thot 
                                        option
  --Thot arg (=0)                       Minimum photosphere of irradiation 
                                        temperature of the outer edge of the 
                                        hot disk, degrees Kelvin. For details 
                                        see --boundcond description
  --F0 arg (=1e+36)                     Initial viscous torque on outer 
                                        boundary of the disk, cgs
  --Mdot0 arg (=0)                      Initial mass accretion rate, g/s. If 
                                        both --F0 and --Mdot0 are specified 
                                        then --Mdot0 is used. Works only when 
                                        --initialcond is setted to sinusF or 
                                        quasistat
  --initialcond arg (=power)            Initial condition viscous torque F or 
                                        surface density Sigma
                                        
                                        Values:
                                          powerF: F ~ xi^powerorder, powerorder
                                        is specified by --powerorder option
                                          powerSigma: Sigma ~ xi^powerorder, 
                                        powerorder is specified by --powerorder
                                        option
                                          sinusF: F ~ sin( xi * pi/2 )
                                          quasistat: F ~ f(xi), where f is 
                                        quasistationary solution found in 
                                        Lipunova, Shakura 2000. f(xi=0) = 0, 
                                        df/dxi(xi=1) = 0
                                        
                                        Here xi is (h - h_in) / (h_out - h_in)
  --powerorder arg (=6)                 Parameter of the powerlaw initial 
                                        condition distributions. This option 
                                        works only with --initialcond=powerF 
                                        and =powerSigma

Parameters of X-ray emission:
  --Cirr arg (=0)                       Irradiation factor
  --irrfactortype arg (=const)          Type of irradiation factor Cirr: const 
                                        (doesn't depend on disk shape, [rad. 
                                        flux] = Cirr  L / [4 pi r^2]), square 
                                        (disk has polynomial shape, [rad. flux]
                                        = Cirr L / [4 pi r^2] [z/r]^2 )
  --dilution arg (=1.7)                 Dilution parameter
  --numin arg (=1)                      Lower bound of X-ray band, keV
  --numax arg (=12)                     Upper bound of X-ray band, keV

Parameters for optical magnitudes calculation:
  --distance arg (=10)                  Distance to the system, kpc

Parameters of disc evolution calculation:
  -T [ --time ] arg (=25)               Computation time, days
  --tau arg (=0.25)                     Time step, days
  --Nx arg (=1000)                      Size of calculation grid
  --gridscale arg (=log)                Type of grid for angular momentum h: 
                                        log or linear
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

License
-------

Copyright (c) 2016, Konstantin L. Malanchev & Galina V. Lipunova.

`Freddi` is distributed under the terms of the GPLv3.

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2016) (in prep.)
