`Freddi` — compute FRED-like light curves of LMXB
=================================================

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. This
code was developed to simulate fast rise exponential decay (FRED) light curves
of low mass X-ray binaries (LMXB) for the paper Lipunova & Malanchev (2016) (in
prep.).

Installation and Usage
----------------------

### Requirements

-   [Boost.Program\_options](http://www.boost.org/doc/libs/release/doc/html/program_options.html)

-   C++ compiler with C++11 support, e.g. `gcc` version 4.8+ or `clang` 3.4+

`Freddi` was tested on Linux and macOS but it should work on Windows as well.

### Installation

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
make install
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Usage

All you need to calculate your first FRED is to type `./freddi` or `freddi` if
it is in your `$PATH`.

`Freddi` has a number of command line options, full list of them with
descriptions can be read using `--help` option. `Freddi` always outputs
`sum.dat` file with temporal dependence of disk values. If `--fulldata` is
specified then `%d.dat` files for each time step are outputted to the same
directory with variuos radial distributions (temperature, viscous torque,
surface density, etc.).

License
-------

Copyright (c) 2016, Konstantin L. Malanchev & Galina V. Lipunova.

`Freddi` is distributed under the terms of the GPLv2 and GPLv3.

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2016) (in prep.)
