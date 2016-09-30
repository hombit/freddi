`Freddi` — compute FRED-like light curves of LMXB
=================================================

Overview
--------

The code solves 1-D evolution equation of Shakura-Sunyaev accretion disk. This
code

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

`Freddi` has a number of command line options, full list with descriptions read
using `--help` option.
