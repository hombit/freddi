`Freddi` — compute FRED-like light curves of LMXB
=================================================

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. This
code was developed to simulate fast rise exponential decay (FRED) light curves
of low mass X-ray binaries (LMXB) for the paper Lipunova & Malanchev (2016) (in
prep.).

Installation
------------

### Docker

If you are Docker user then you can run `Freddi` as Docker executable, skip all
installation instructions and go straight to the Usage section

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
`--dir` option.

See full list of command line options with `--help` option.

License
-------

Copyright (c) 2016, Konstantin L. Malanchev & Galina V. Lipunova.

`Freddi` is distributed under the terms of the GPLv3.

Please, accompany any results obtained using this code with reference to
Lipunova & Malanchev (2016) (in prep.)
