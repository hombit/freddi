from glob import glob

from setuptools import setup, Extension
from Cython.Build import cythonize

cpp_source = glob('src/*.cpp')
cpp_source.remove('src/main.cpp')

extensions = [
    Extension('freddi.args',
              cpp_source + ['freddi/args.pyx'],
              extra_compile_args=['-std=c++11'],
              include_dirs=['src/'],
              libraries=['boost_program_options'])
]

setup(ext_modules=cythonize(extensions))
