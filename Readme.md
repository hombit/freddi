# `Freddi` — compute FRED-like light curves of LMXB

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Physical Background](#physical-background) 
- [Accretion disk wind](#accretion-disk-wind) 
- [Development guide](#development-guide)
- [Questions and comments](#questions-and-comments)
- [License](#license)
- [BibTex](#bibtex)

## Overview

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. The
code is developed to simulate fast-rise exponential-decay (FRED) light curves of
low mass X-ray binaries (LMXBs) for the paper “Determination of the turbulent
parameter in the accretion disks: effects of self-irradiation in 4U 1543-47
during the 2002 outburst” by Lipunova & Malanchev (2017)
[2017MNRAS.468.4735L](http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L).

`Freddi` is written on C++ and available as a couple of binary executables and
a Python module.

Note that the original `Freddi` version 1 introduced in Lipunova & Malanchev (2017)
[2017MNRAS.468.4735L](http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L) is still
available in the [`v1` git branch](https://github.com/hombit/freddi/tree/v1).

## Installation

### Executables

`Freddi` is represented by two binary executables: the black hole version
`freddi` and the neutron star version `freddi-ns`.

#### Docker

If you are familiar with [Docker](http://docker.com) then you can use
pre-compiled binaries inside Docker container:
```sh
docker run -v "`pwd`":/data --rm -ti ghcr.io/hombit/freddi freddi -d/data
docker run -v "`pwd`":/data --rm -ti ghcr.io/hombit/freddi freddi-ns --Bx=1e8 -d/data
```

#### Build from source

`Freddi` has following build dependencies:

- [Boost](http://www.boost.org/) 1.57+
- [CMake](https://cmake.org) with a back-end build system like Make or Ninja
- C++ compiler with C++17 support, e.g. `gcc` version 8+ or `clang` 5+

Get requirements on Debian based systems (e.g. Ubuntu):
```sh
apt-get install g++ cmake libboost-all-dev
```

On Red-Hat based systems (e.g. Fedora):
```sh
dnf install gcc-c++ cmake boost-devel
```

On macOS via [Homebrew](https://brew.sh):
```sh
brew install cmake boost
```

Get `Freddi` source code and compile it:
```sh
git clone https://github.com/hombit/freddi
cd freddi
mkdir cmake-build
cd cmake-build
cmake .. # -DSTATIC_LINKING=TRUE
cmake --build .
```

Uncomment `-DSTATIC_LINKING=TRUE` to link against static Boost libraries

Now you should have both `freddi` and `freddi-ns` executables in the build
directory. You can install these binaries and the default configuration
file `freddi.ini` by running
```sh
cmake --install . --prefix=PREFIX  # replace with preferable location
```

`Freddi` is known to be built on Linux and macOS.

### Python

[![PyPI version](https://badge.fury.io/py/freddi.svg)](https://badge.fury.io/py/freddi)

Python 2 isn't supported, use Python 3 instead.

`Freddi` pre-compiled x86-64 Linux packages for several Python versions
are available on <https://pypi.org/project/freddi/> and can be used as is,
while for other configurations you should have C++ compiler and Boost
libraries in your system before running this command:

```sh
# Please upgrade your pip
python3 -m pip install -U pip
# Depending on your Python setup, you need or need not --user flag
python3 -m pip install --user freddi
```

[`astropy`](https://astropy.org) is an optional requirement which must be
installed to use dimensional input via `Freddi.from_astropy`

## Usage

### Executables

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

#### <a name="usage-executables-options"> Options

The full list of command line options is viewed with `--help` option. Default
values are given in brackets.

```sh
./freddi --help
```

<details><summary>expand</summary>

```
Freddi: numerical calculation of accretion disk evolution:

General options:
  -h [ --help ]                    Produce help message
  --config arg                     Set additional configuration filepath
  --prefix arg (=freddi)           Set prefix for output filenames. Output file
                                   with distribution of parameters over time is
                                   PREFIX.dat
  --stdout                         Output temporal distribution to stdout 
                                   instead of PREFIX.dat file
  -d [ --dir ] arg (=.)            Choose the directory to write output files. 
                                   It should exist
  --precision arg (=12)            Number of digits to print into output files
  --tempsparsity arg (=1)          Output every k-th time moment
  --fulldata                       Output files PREFIX_%d.dat with radial 
                                   structure for every time step. Default is to
                                   output only PREFIX.dat with global disk 
                                   parameters for every time step

Basic binary and disk parameter:
  -a [ --alpha ] arg               Alpha parameter of Shakura-Sunyaev model
  --alphacold arg                  Alpha parameter of cold disk, currently it 
                                   is used only for Sigma_minus, see 
                                   --Qirr2Qvishot. Default is --alpha values 
                                   divided by ten
  -M [ --Mx ] arg                  Mass of the central object, in the units of 
                                   solar masses
  --kerr arg (=0)                  Dimensionless Kerr parameter of the black 
                                   hole
  --Mopt arg                       Mass of the optical star, in units of solar 
                                   masses
  --rochelobefill arg (=1)         Dimensionless factor describing a size of 
                                   the optical star. Polar radius of the star 
                                   is rochelobefill * (polar radius of critical
                                   Roche lobe)
  --Topt arg (=0)                  Thermal temperature of the optical star, in 
                                   units of kelvins
  -P [ --period ] arg              Orbital period of the binary system, in 
                                   units of days
  --rin arg                        Inner radius of the disk, in the units of 
                                   the gravitational radius of the central 
                                   object GM/c^2. If it isn't set then the 
                                   radius of ISCO orbit is used defined by --Mx
                                   and --kerr values
  -R [ --rout ] arg                Outer radius of the disk, in units of solar 
                                   radius. If it isn't set then the tidal 
                                   radius is used defined by --Mx, --Mopt and 
                                   --period values
  --risco arg                      Innermost stable circular orbit, in units of
                                   gravitational radius of the central object 
                                   GM/c^2. If it isn't set then the radius of 
                                   ISCO orbit is used defined by --Mx and 
                                   --kerr values

Parameters of the disk mode:
  -O [ --opacity ] arg (=Kramers)  Opacity law: Kramers (varkappa ~ rho / 
                                   T^7/2) or OPAL (varkappa ~ rho / T^5/2)
  --Mdotout arg (=0)               Accretion rate onto the disk through its 
                                   outer radius
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
  --Qirr2Qvishot arg (=0)          Minimum Qirr / Qvis ratio at the outer edge 
                                   of the hot disk to switch evolution from 
                                   temperature-based regime to 
                                   Sigma_minus-based regime (see Eq. A.1 in 
                                   Lasota et al. 2008, --alphacold value is 
                                   used as alpha parameter)
  --initialcond arg (=powerF)      Type of the initial condition for viscous 
                                   torque F or surface density Sigma
                                   
                                   Values:
                                     [xi = (h - h_in) / (h_out - h_in)]
                                     powerF: F ~ xi^powerorder, powerorder is 
                                   specified by --powerorder option
                                     linearF: F ~ xi, specific case of powerF 
                                   but can be normalised by --Mdot0, see its 
                                   description for details  powerSigma: Sigma ~
                                   xi^powerorder, powerorder is specified by 
                                   --powerorder option
                                     sineF: F ~ sin( xi * pi/2 )
                                     gaussF: F ~ exp(-(xi-mu)**2 / 2 sigma**2),
                                   mu and sigma are specified by --gaussmu and 
                                   --gausssigma options
                                     quasistat: F ~ f(h/h_out) * xi * h_out/h, 
                                   where f is quasi-stationary solution found 
                                   in Lipunova & Shakura 2000. f(xi=0) = 0, 
                                   df/dxi(xi=1) = 0
                                   
  --F0 arg                         Initial maximum viscous torque in the disk, 
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
                                   linearF, sinusF or quasistat
  --powerorder arg                 Parameter for the powerlaw initial condition
                                   distribution. This option works only with 
                                   --initialcond=powerF or powerSigma
  --gaussmu arg                    Position of the maximum for Gauss 
                                   distribution, positive number not greater 
                                   than unity. This option works only with 
                                   --initialcond=gaussF
  --gausssigma arg                 Width of for Gauss distribution. This option
                                   works only with --initialcond=gaussF
  --windtype arg (=no)             Type of the wind
                                   
                                     no: no wind
                                     SS73C: super-Eddington spherical wind from
                                   Shakura-Sunyaev 1973
                                     ShieldsOscil1986: toy wind model from 
                                   Shields et al. 1986 which was used to obtain
                                   oscillations in the disk luminosity. 
                                   Requires --windC_w and --windR_w to be 
                                   specified
                                     Janiuk2015: super-Eddington wind from 
                                   Janiuk et al 2015. Requires --windA_0 and 
                                   --windB_1 to be specified
                                     Shields1986: thermal wind from Begelman et
                                   al. 1983 and Shields et al. 1986. Requires 
                                   --windXi_max, --windT_ic and --windPow to be
                                   specified
                                     Woods1996AGN: thermal AGN wind from Woods 
                                   et al. 1996. Requires --windC_0 and 
                                   --windT_ic to be specified
                                     Woods1996: thermal wind from Woods et al. 
                                   1996. Requires --windXi_max, --windT_ic and 
                                   --windPow to be specified
                                     toy: a toy wind model used in 
                                   arXiv:2105.11974, the mass loss rate is 
                                   proportional to the central accretion rate. 
                                   Requires --windC_w to be specified
                                   
  --windC_w arg                    The ratio of the mass loss rate due to wind 
                                   to the central accretion rate, |Mwind|/Macc
  --windR_w arg                    The ratio of the wind launch radius to the 
                                   outer disk radius, Rwind/Rout
  --windA_0 arg                    Dimensionless parameter characterizing the 
                                   strength of the super-Eddington wind in the 
                                   framework of the model Janiuk et al. 2015. 
                                   Effective value range from 10 to 25
  --windB_1 arg                    The quantity is of the order of unity. 
                                   Characterizes the relationship between the 
                                   change in energy per particle and virial 
                                   energy.
                                   E = B_1 * k * T
  --windXi_max arg                 Ionization parameter, the ratio of the 
                                   radiation and gas pressures
  --windT_ic arg                   Inverse Compton temperature, K. 
                                   Characterizes the hardness of the 
                                   irradiating spectrum
  --windPow arg                    Multiplicative coefficient to control wind 
                                   power
  --windC_0 arg                    Characteristic column density of the wind 
                                   mass loss rate from Woods et al. 1996 model,
                                   g/(s*cm^2). For AGN approx value is 3e-13 
                                   g/(s*cm^2)

Parameters of self-irradiation.
Qirr = Cirr * (H/r / 0.05)^irrindex * L * psi / (4 pi R^2), where psi is angular distrbution of X-ray radiation:
  --Cirr arg (=0)                  Irradiation factor for the hot disk
  --irrindex arg (=0)              Irradiation index for the hot disk
  --Cirrcold arg (=0)              Irradiation factor for the cold disk
  --irrindexcold arg (=0)          Irradiation index for the cold disk
  --h2rcold arg (=0)               Seme-height to radius ratio for the cold 
                                   disk, it affects disk shadow in star
  --angulardistdisk arg (=plane)   Angular distribution of the disk X-ray 
                                   radiation. Values: isotropic, plane

Parameters of flux calculation:
  --colourfactor arg (=1.7)        Colour factor to calculate X-ray flux
  --emin arg (=1)                  Minimum energy of X-ray band, keV
  --emax arg (=12)                 Maximum energy of X-ray band, keV
  --staralbedo arg (=0)            Part of X-ray radiation reflected by optical
                                   star, (1 - albedo) heats star's photosphere.
                                   Used only when --starflux is specified
  -i [ --inclination ] arg (=0)    Inclination of the system, degrees
  --ephemerist0 arg (=0)           Ephemeris for the time of the minimum of the
                                   orbital light curve T0, phase zero 
                                   corresponds to inferior conjunction of the 
                                   optical star, days
  --distance arg                   Distance to the system, kpc
  --colddiskflux                   Add Fnu for cold disk into output file. 
                                   Default output is for hot disk only
  --starflux                       Add Fnu for irradiated optical star into 
                                   output file. See --Topt, --starlod and 
                                   --h2rcold options. Default is output for the
                                   hot disk only
  --lambda arg                     Wavelength to calculate Fnu, Angstrom. You 
                                   can use this option multiple times. For each
                                   lambda one additional column with values of 
                                   spectral flux density Fnu [erg/s/cm^2/Hz] is
                                   produced
  --passband arg                   Path of a file containing tabulated 
                                   passband, the first column for wavelength in
                                   Angstrom, the second column for transmission
                                   factor, columns should be separated by 
                                   spaces

Parameters of disk evolution calculation:
  --inittime arg (=0)              Initial time moment, days
  -T [ --time ] arg                Time interval to calculate evolution, days
  --tau arg                        Time step, days
  --Nx arg (=1000)                 Size of calculation grid
  --gridscale arg (=log)           Type of grid for angular momentum h: log or 
                                   linear
  --starlod arg (=3)               Level of detail of the optical star 3-D 
                                   model. The optical star is represented by a 
                                   triangular tile, the number of tiles is 20 *
                                   4^starlod


```
</details>

```sh
./freddi-ns --help
```

<details><summary>expand</summary>

```
Freddi NS: numerical calculation of accretion disk evolution:

General options:
  -h [ --help ]                         Produce help message
  --config arg                          Set additional configuration filepath
  --prefix arg (=freddi)                Set prefix for output filenames. Output
                                        file with distribution of parameters 
                                        over time is PREFIX.dat
  --stdout                              Output temporal distribution to stdout 
                                        instead of PREFIX.dat file
  -d [ --dir ] arg (=.)                 Choose the directory to write output 
                                        files. It should exist
  --precision arg (=12)                 Number of digits to print into output 
                                        files
  --tempsparsity arg (=1)               Output every k-th time moment
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
                                          [xi = (h - h_in) / (h_out - h_in)]
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
                                          quasistat-ns: ???
                                        
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
  --windtype arg (=no)                  Type of the wind
                                        
                                          no: no wind
                                          SS73C: super-Eddington spherical wind
                                        from Shakura-Sunyaev 1973
                                          ShieldsOscil1986: toy wind model from
                                        Shields et al. 1986 which was used to 
                                        obtain oscillations in the disk 
                                        luminosity. Requires --windC_w and 
                                        --windR_w to be specified
                                          Janiuk2015: super-Eddington wind from
                                        Janiuk et al 2015. Requires --windA_0 
                                        and --windB_1 to be specified
                                          Shields1986: thermal wind from 
                                        Begelman et al. 1983 and Shields et al.
                                        1986. Requires --windXi_max, --windT_ic
                                        and --windPow to be specified
                                          Woods1996AGN: thermal AGN wind from 
                                        Woods et al. 1996. Requires --windC_0 
                                        and --windT_ic to be specified
                                          Woods1996: thermal wind from Woods et
                                        al. 1996. Requires --windXi_max, 
                                        --windT_ic and --windPow to be 
                                        specified
                                          toy: a toy wind model used in 
                                        arXiv:2105.11974, the mass loss rate is
                                        proportional to the central accretion 
                                        rate. Requires --windC_w to be 
                                        specified
                                        
  --windC_w arg                         The ratio of the mass loss rate due to 
                                        wind to the central accretion rate, 
                                        |Mwind|/Macc
  --windR_w arg                         The ratio of the wind launch radius to 
                                        the outer disk radius, Rwind/Rout
  --windA_0 arg                         Dimensionless parameter characterizing 
                                        the strength of the super-Eddington 
                                        wind in the framework of the model 
                                        Janiuk et al. 2015. Effective value 
                                        range from 10 to 25
  --windB_1 arg                         The quantity is of the order of unity. 
                                        Characterizes the relationship between 
                                        the change in energy per particle and 
                                        virial energy.
                                        E = B_1 * k * T
  --windXi_max arg                      Ionization parameter, the ratio of the 
                                        radiation and gas pressures
  --windT_ic arg                        Inverse Compton temperature, K. 
                                        Characterizes the hardness of the 
                                        irradiating spectrum
  --windPow arg                         Multiplicative coefficient to control 
                                        wind power
  --windC_0 arg                         Characteristic column density of the 
                                        wind mass loss rate from Woods et al. 
                                        1996 model, g/(s*cm^2). For AGN approx 
                                        value is 3e-13 g/(s*cm^2)

Parameters of accreting neutron star:
  --nsprop arg (=dummy)                 Neutron star geometry and radiation 
                                        properties name and specifies default 
                                        values of --Rx, --Risco and --freqx
                                        
                                        Values:
                                          dummy: NS radiation efficiency is R_g
                                        * (1 / R_x - 1 / 2R_in), default 
                                        --freqx is 0, default Rx is 1e6, 
                                        default Risco is Kerr value
                                          newt: NS radiation efficiency is a 
                                        functions of NS frequency, that's why 
                                        --freqx must be specified explicitly
                                          sibgatullin-sunyaev2000: NS radiation
                                        efficiency, and default values of Rx 
                                        and Risco are functions of NS 
                                        frequency, that's why --freqx must be 
                                        specified explicitly
  --freqx arg                           Accretor rotation frequency, Hz. This 
                                        parameter is not linked to --kerr, 
                                        agree them yourself
  --Rx arg                              Accretor radius, cm
  --Bx arg                              Accretor polar magnetic induction, G
  --hotspotarea arg (=1)                ??? Relative area of hot spot on the 
                                        accretor
  --epsilonAlfven arg (=1)              Factor in Alfven radius formula
  --inversebeta arg (=0)                ???
  --Rdead arg (=0)                      Maximum inner radius of the disk that 
                                        can be obtained, it characterises 
                                        minimum torque in the dead disk, cm
  --fptype arg (=no-outflow)            ??? Accretor Mdot fraction mode fp.
                                        
                                        Values:
                                          no-outflow: all the matter passing 
                                        inner disk radius falling onto neutron 
                                        star, fp = 1
                                          propeller: all the matter flows away 
                                        from both disk and neutron star, fp = 0
                                          corotation-block: like 'no-otflow' 
                                        when Alfven radius is smaller than 
                                        corotation radius, like 'propeller' 
                                        otherwise
                                          geometrical: generalisation of 
                                        'corotation-block' for the case of not 
                                        co-directional of disk rotation axis 
                                        and neutron star magnetic field axis. 
                                        Requires --fp-geometrical-chi to be 
                                        specified
                                          eksi-kutlu2010: ???
                                          romanova2018: ???, requires 
                                        --romanova2018-par1 and 
                                        --romanova2018-par2 to be specified
  --fp-geometrical-chi arg              angle between disk rotation axis and 
                                        neutron star magnetic axis for 
                                        --fptype=geometrical, degrees
  --romanova2018-par1 arg               ??? par1 value for --fptype=romanova201
                                        8 and --kappattype=romanova2018
  --romanova2018-par2 arg               ??? par2 value for --fptype=romanova201
                                        8 and --kappattype=romanova2018
  --kappattype arg (=const)             kappa_t describes how strong is 
                                        interaction between neutron star 
                                        magnitosphere and disk, magnetic torque
                                        is kappa_t(R) * mu^2 / R^3. This 
                                        parameter describes type of radial 
                                        destribution of this parameter
                                        
                                        Values:
                                          const: doesn't depend on radius, 
                                        kappa_t = value. Requires 
                                        --kappat-const-value to be specified
                                          corstep: kappa_t is 'in' inside 
                                        corotation radius, and 'out' outside. 
                                        Requires --kappat-corstep-in and 
                                        --kappat-corstep-out to be specified
                                          romanova2018: similar to corstep 
                                        option, but the outside value is 
                                        reduced by the portion taken away by 
                                        wind (see Table 2 of 
                                        Romanova+2018,NewA,62,94). Requires 
                                        --kappat-romanova2018-in, 
                                        --kappat-romanova2018-out 
                                        --romanova2018-par1 and --romanova-par2
                                        to be specified
  --kappat-const-value arg (=0.33333333333333331)
                                        kappa_t value for --kappattype=const
  --kappat-corstep-in arg (=0.33333333333333331)
                                        kappa_t value inside corotation radius 
                                        for --kappattype=corstep
  --kappat-corstep-out arg (=0.33333333333333331)
                                        kappa_t value outside corotation radius
                                        for --kappattype=corstep
  --kappat-romanova2018-in arg (=0.33333333333333331)
                                        kappa_t value inside corotation radius 
                                        for --kappattype=romanova2018
  --kappat-romanova2018-out arg (=0.33333333333333331)
                                        kappa_t value outside corotation radius
                                        for --kappattype=romanova2018
  --nsgravredshift arg (=off)           Neutron star gravitational redshift 
                                        type.
                                        
                                        Values:
                                          off: gravitational redshift doesn't 
                                        taken into account
                                          on: redshift is (1 - R_sch / Rx), 
                                        where R_sch = 2GM/c^2

Parameters of self-irradiation.
Qirr = Cirr * (H/r / 0.05)^irrindex * L * psi / (4 pi R^2), where psi is angular distrbution of X-ray radiation:
  --Cirr arg (=0)                       Irradiation factor for the hot disk
  --irrindex arg (=0)                   Irradiation index for the hot disk
  --Cirrcold arg (=0)                   Irradiation factor for the cold disk
  --irrindexcold arg (=0)               Irradiation index for the cold disk
  --h2rcold arg (=0)                    Seme-height to radius ratio for the 
                                        cold disk, it affects disk shadow in 
                                        star
  --angulardistdisk arg (=plane)        Angular distribution of the disk X-ray 
                                        radiation. Values: isotropic, plane
  --angulardistns arg (=isotropic)      Angular distribution type of the 
                                        neutron star X-ray radiation. Values: 
                                        isotropic, plane

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
  --starflux                            Add Fnu for irradiated optical star 
                                        into output file. See --Topt, --starlod
                                        and --h2rcold options. Default is 
                                        output for the hot disk only
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


```
</details>

**Write which options are mandatory**

Also you can use `freddi.ini` configuration file to store options. This [INI
file](https://en.wikipedia.org/wiki/INI_file) contains lines `option=value`,
option names are the as provided by the help message above. Command line option
overwrites configuration file option. For example, [see
default](https://github.com/hombit/freddi/blob/master/freddi.ini) `freddi.ini`.

Paths where this file is searched are `./freddi.ini` (execution path),
`$HOME/.config/freddi/freddi.ini`, `/usr/local/etc/freddi.ini` and
`/etc/freddi.ini`. You can provide configuration file to Docker container as a
volume: `` -v "`pwd`/freddi.ini":/etc/freddi.ini ``.

#### Output values

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

#### <a name="usage-executables-example"></a> Example

The following arguments instruct `Freddi` to calculate the decay of the outburst
in the disk with the constant outer radius equal to 1 solar radius. The Kerr
black hole at the distance of 5 kpc has the mass of 9 solar masses, and the
Kerr's parameter is 0.4. The outer disk is irradiated with Cirr=1e-3.
**Discuss all options used in the example***

```sh
./freddi --alpha=0.5 --Mx=9 --rout=1 --period=0.5 --Mopt=0.5 --time=50 \
  --tau=0.25 --dir=data/ --F0=2e+37 --colourfactor=1.7 --Nx=1000 \
  --distance=5 --gridscale=log --kerr=0.4 --Cirr=0.001 --opacity=OPAL \
  --initialcond=quasistat --windtype=Woods1996 --windXi_max=10 --windT_ic=1e8 \
  --windPow=1 
```

### Python

Python bindings can be used as a convenient way to run and analyse Freddi
simulations.

#### Initializing

You can prepare simulation set-up initializing `Freddi`
(for black hole accretion disk) or `FreddiNeutronStar` (for NS) class instance.
These classes accept keyword-only arguments which have the same names and
meanings as [command line options](#usage-executables-options), but with three
major exceptions:
 1. Python package doesn't provide any file output functionality, that's why output arguments like `config`, `dir`, `fulldata`, `starflux`, `lambda` or `passband` are missed;
 2. all values are assumed to be in CGS units, but you can use `Freddi.from_asrtopy` for dimensional values (see details bellow);
 3. parameters of wind, NS `fp` and NS `kappa` models are passed as dictionaries (see specification bellow).

The following code snippet would set-up roughly the same simulation as
[the command-line example](#usage-executables-example)

```python
from freddi import Freddi

freddi = Freddi(
    alpha=0.5, Mx=9*2e33, rout=1*7e10, period=0.5*86400, Mopt=0.5*2e33,
    time=50*86400, tau=0.25*86400, F0=2e+37, colourfactor=1.7, Nx=1000,
    distance=5*3e21, gridscale='log', kerr=0.4, Cirr=0.001, opacity='OPAL',
    initialcond='quasistat', wind='Woods1996',
    windparams=dict(Xi_max=10, T_iC=1e8, W_pow=1),
)
```

Alternatively we can do the same using `from_astropy` class-method which casts
all [`astropy.units.Quantity`](https://docs.astropy.org/en/stable/units/quantity.html)
objects to CGS values. Note that dimensionality isn't checked, and technically
it just does `arg.cgs.value` for every `Quantity` argument.

```python
import astropy.units as u
from freddi import Freddi

freddi = Freddi.from_astropy(
    alpha=0.5, Mx=9*u.Msun, rout=1*u.Rsun, period=0.5*u.day, Mopt=0.5*u.Msun,
    time=50*u.day, tau=0.25*u.day, F0=2e+37, colourfactor=1.7, Nx=1000,
    distance=5*u.kpc, gridscale='log', kerr=0.4, Cirr=0.001, opacity='OPAL',
    initialcond='quasistat', wind='Woods1996',
    windparams=dict(Xi_max=10, T_iC=1e8*u.K, W_pow=1),
)
```

Wind model parameters are specified by `windparams` argument which should be
a `dict` instance with string keys and numeric values. Command option to 
`windparams` keys relation is: `--windC_w -> C_w`, `--windR_w -> R_w`,
`--windA_0 -> A_0`, `--windB_1 -> B_1`, `--windXi_max -> Xi_max`,
`windT_ic -> T_ic`, `--windPow -> Pow`, `windC_0 -> C_0`.

Neutron star f_p model parameters are specified by `fpparams` mapping with the
same structure as `windparams`. Command options to `fpparams` keys relation is:
`--fp-geometrical-chi -> chi`, `romanova2018-par1 -> par1`,
`romanova2018-par2 -> par2`.

Neutron star kappa_t model parameters are specified by `kappatparams` mapping
with the same structure as `windparams`. Command options to the mapping keys
relation is: `--kappat-const-value -> value`, `--kappat-corstep-in -> in`,
`kappat-corstep-out -> out`, `--kappat-romanova2018-in -> in`,
`--kappat-romanova2018-out -> out`, `romanova2018-par1 -> par1`,
`--romanova2018-par2 -> par2`

#### Running

There are two ways to run a simulation: iterating over time steps, and run the
whole simulation in one shot. Note that in both cases your `Freddi` object is
mutating and represents the current state of the accretion disk.

Here we use iterator interface which yields another `Freddi` object for each
time moment.

```python
import astropy.units as u
from freddi import Freddi

freddi = Freddi.from_astropy(
    alpha=0.5, Mx=9*u.Msun, rout=1*u.Rsun, period=0.5*u.day, Mopt=0.5*u.Msun,
    time=20*u.day, tau=1.0*u.day, Mdot0=5e18, distance=10*u.kpc,
    initialcond='quasistat',
)

for state in freddi:
    print(f't = {state.t} s, Mdot = {state.Mdot:g} g/s')
assert state.t == freddi.t
```

In this example we run a simulation via `.evolve()` method which returns
`EvolutionResult` object keeping all evolution states internally and providing
temporal distribution of disk's properties.

```python
import astropy.units as u
import matplotlib.pyplot as plt
from freddi import Freddi

freddi = Freddi.from_astropy(
    alpha=0.5, Mx=9*u.Msun, rout=1*u.Rsun, period=0.5*u.day, Mopt=0.5*u.Msun,
    time=20*u.day, tau=1.0*u.day, Mdot0=5e18, distance=10*u.kpc,
    initialcond='quasistat',
)

result = freddi.evolve()
assert result.t[-1] == freddi.t

# Plot Mdot(t)
plt.figure()
plt.title('Freddi disk evolution: accretion rate')
plt.xlabel('t, day')
plt.ylabel(r'$\dot{M}$, g/cm')
plt.plot(result.t / 86400, result.Mdot)
plt.show()

# Plot all F(h) profiles
plt.figure()
plt.title('Freddi disk evolution: viscous torque')
plt.xlabel('r, cm')
plt.ylabel('F, dyn cm')
plt.xscale('log')
plt.yscale('log')
plt.plot(result.R.T, result.F.T)
plt.show()

# Plot evolution of effective temperature of the outer hot disk ring
plt.figure()
plt.title('Freddi disk evolution: outer effective temperature')
plt.xlabel('t, day')
plt.ylabel('T, K')
plt.plot(result.t / 86400, result.last_Tph)
plt.show()
```

#### Properties and methods

`Freddi`, `FreddiNeutronStar` and `EvolutionResult` objects contain dozens of
properties returning various physical values like `t` for time moment,
`Mdot` for accretion rate onto central object, `R` for radius, `F` for torque,
`Tph` for effective temperature and so on. `first_*` and `last_*` properties
are used to access innermost and outermost values of radial-distributed
quantities. The complete list of properties can be obtained by `dir(Freddi)` or
`dir(FreddiNeutronStar)`. Note that the most properties are lazy-evaluated and
require some time to access first time. `EvolutionRadius` provides all the
same properties as underlying `Freddi` or `FreddiNeutronStar` objects but with
additional array dimension for temporal distribution, so if `Freddi.Lx` is a
scalar then `EvolutionResult.Lx` is 1-D `numpy` array of `(Nt,)` shape,
if `Freddi.Sigma` is 1-D array of `(Nx,)` shape, then
`EvolutionResult.Sigma` is 2-D array of `(Nt, Nx)` shape. Also, note that if
disk shrinks during a simulation, the missing values of `EvolutionResult`
properties are filled by NaN.

All three classes have `flux(lmbd, region='hot', phase=None)` method which can
be used to find spectral flux density per unit frequency for optical
emission. `lmbd` argument can be a scalar or a multidimensional `numpy` array
of required wavelengths in cm; `region` could be one of "hot" (hot disk),
"cold" (cold disk), "disk" ("hot" + "cold"), "star" (companion star), and
"all" ("hot" + "cold" + "star"); `phase` is a binary system orbital phase in
radians, it is required for `region="star"` and `region="all"` only, it can be
calculated as `2π t / period + constant`.

All properties and methods return values in CGS units.

## Physical Background

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


### Accretion disk wind

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
`--windtype` option. The thermal wind model ([Woods et al. 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...461..767W)),
which implies that the outflow of matter occurs due to the heating of the outer parts of the disk
by a central radiation source, can be chosen by setting `--windtype=Woods1996`. 
The option `--windtype=Janiuk15` corresponds to the model from work [Janiuk et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..92J)
where the wind is started in the super-Eddington regime.
When choosing option `--windtype=Janiuk15`, the you must also specify the values of
the super-Eddington wind parameters with `--windA0` and `--windB1` options.
You can also select the `--windtype=toy` option, which corresponds to a toy wind model when the user sets 
the wind strength relatively to the accretion rate using the option `--windPow`.

#### Compton-heated wind

At the moment, `Freddi` is more focused on simulating outbursts taking into account the thermal wind (`--windtype=Woods1996` option). 
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

Choosing option `--windtype=Woods1996`, it is necessary to set the value of the ionization parameter Xi
(which is proportional to the ratio of the radiation and gas pressures) by the option `--windXi_max` and the Compoton temperature T_ic 
(which determines the hardness of the irradiating spectrum and the size of the region where the wind operates) by the option `--windT_ic`. 

### Companion star irradiation

We use a simple model of irradiated star to simulate periodic variability and
X-ray thermalization by a companion's photosphere. Our model assumes that the
companion star's shape corresponds to equipotential surface which size is set
by `--rochelobefill` option, unity means that star fills its Roche lobe, any
smaller value decreases star's polar radius correspondingly. Technically,
star's surface is built from `20 * 4^starlod` triangles, use `--starlod` to
set level of detail, `--starlod=3` should give few percent precision. Every
triangle has black-body spectrum with bolometric luminosity given by a sum of
star's own luminosity (set by `--Topt`) and irradiation flux multiplied by 
unity minus albedo (set by `--staralbedo`).

Please note that the model is limited and doesn't implement limb darking or
eclipsing.

## Development guide

### Source code and tests

`Freddi` uses [Cmake](https://cmake.org) as a build system.

The C++ source code is located in `cpp` folder which has following structure:
- `main.cpp` and `main-ns.cpp` implements `main()` function for `freddi` and `freddi-ns` correspondingly;
- `include` for library header files, it has `ns` sub-folder for neutron star related stuff;
- `src` for library C++ files, it also has `ns` sub-folder;
- `test` provides library unit tests;
- `pywrap` has both header and source files for `Boost::Python`/`Boost::NumPy` bindings.

Note, that we require C++17 standard (while not having idiomatic C++17 code),
and require code to be compiled by modern GCC and CLang on Linux. Please write
unit tests where you can and use `ctest` to check they pass.

The Python project is specified by `pyproject.toml` (which just lists build
requirements), `setup.py` and `MANIFEST.in` files, we use
[`scikit-build`](https://scikit-build.readthedocs.io/) as a build system.
`scikit-build` uses Python-related section of `CMakeLists.txt` to build C++
source code into Python extension, and accomplish it with Python files located
in `python/freddi` directory. Use `python setup.py build_ext` to build the
extension, optionally with `-DSTATIC_LINKING=TRUE` to link `Boost::Filesystem`,
`Boost::Python` and `Boost::NumPy` statically. Please, pay attention to two
last libraries, because they should be built against the same Python version as
you use.

`python/test` contains some tests, you can run them by `python3 setup.py test`.
- `test_freddi.py` and `test_ns.py` contain unit tests for Python source;
- `test_analytical.py` contains integration tests to compare analytical solutions of the equation of disk viscous evolution with the numerical solutions of `Freddi`;
- `regression.py` contains regression tests to be sure that 1) the `Freddi` output is stable, and 2) the Python code gives the same results as binary executables do.

The regression test data are located in `python/test/data`. Sometimes you need
to update these regression data, for example when you introduce new
command-line option with a default value, add new output column or fix some bug
in physical model. For these purposes you can use `generate_test_data.sh`
script located in this folder.

`Dockerfile` is used to build a Docker image with statically-linked binaries,
and `Dockerfile.python` is used to build a Docker image with
[`manylinux`](https://github.com/pypa/manylinux)-compatible Python wheels.


### Continuous integration

We use [Github Actions](https://github.com/hombit/freddi/actions) as a
continuous integration (CI) system. The workflow file is located in
`.github/workflows/main.yml` and a couple of auxiliary files are located in
`.ci` folder. CI allows us to test new commits to prevent different bugs:
- `gcc` and `clang` actions test binaries building, execute sample `Freddi` programs, run C++ unit tests, perform C++ regression tests, and check the consistency of the `Readme.md` with programs' `--help` output
- `cpython` action builds Python extension module and runs all Python tests
- `docker-exe` builds a Docker image using `Dockerfile` and execute sample `Freddi` programs inside a Docker container
- `docker-python` builds a Docker image using `Dockefile.python`, uses wheels it has built to build Python Docker images for several Python versions using `.ci/Dockerfile-test-wheels`, and runs sample Python scripts with `freddi.Freddi` class

### This Readme

Please keep Readme updated. You can update the help messages in the
[Usage](#usage) section using `.ci/update-help-readme.py` script.

### Release new version

Check-list:

- Create `git` tag
- Build new `freddi` image using `Dockerfile`
- Build new `freddi-python` image using `Dockerfile.python`
- Run `docker run --rm -ti docker-python:VERSION sh -c "python3.7 -m twine upload /dist/*"` to upload source code distribution and x86-64 Python wheels onto PyPi.org
- [Optional] Build executables for GitHub release
- [Optional] Build and upload macOS wheels
- [Optional] Build and upload Linux AArch64 wheels
- Crate new GitHub release


## Questions and comments

If you have any problems, questions, or comments, please address them to
[Issues](https://github.com/hombit/freddi/issues) or to hombit\@gmail.com

## License

Copyright (c) 2016–2021, Konstantin L. Malanchev, Galina V. Lipunova & Artur L. Avakyan.

`Freddi` is distributed under the terms of the
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2017) [2017MNRAS.468.4735L](http://adsabs.harvard.edu/abs/2017MNRAS.468.4735L), and for the case of windy calculations please also refer Avakyan et al. (2021) [2021arXiv210511974A](https://ui.adsabs.harvard.edu/abs/2021arXiv210511974A/).

## BibTex
```bibtex
@ARTICLE{2017MNRAS.468.4735L,
   author = { {Lipunova}, G.~V. and {Malanchev}, K.~L.},
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

@ARTICLE{2021arXiv210511974A,
       author = { {Avakyan}, A.~L. and {Lipunova}, G.~V. and {Malanchev}, K.~L. and {Shakura}, N.~I.},
        title = "{Change in the orbital period of a binary system due to an outburst in a windy accretion disc}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena},
         year = 2021,
        month = may,
          eid = {arXiv:2105.11974},
        pages = {arXiv:2105.11974},
archivePrefix = {arXiv},
       eprint = {2105.11974},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv210511974A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
