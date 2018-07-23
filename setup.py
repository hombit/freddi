#!/usr/bin/env python

from glob import glob
from sys import argv

import numpy as np
from setuptools import setup, Extension


try:
    argv.remove('--cython')
    USE_CYTHON = True
except ValueError:
    USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.cpp'
cython_source = 'cython/freddi' + ext

cpp_source = glob('cpp/*.cpp')

extensions = [
    Extension('freddi',
              cpp_source + [cython_source],
              extra_compile_args=['-std=c++11'],
              include_dirs=['cpp', np.get_include()],
              libraries=['boost_program_options'])
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, annotate=True, force=True)


setup(
    name='freddi',
    version='2.0a0.dev0',
    url='http://xray.sai.msu.ru/~malanchev/freddi/',
    license='MIT',
    author='Konstantin Malanchev',
    author_email='malanchev@sai.msu.ru',
    description='Compute FRED light curves of LMXBs outbursts',
    ext_modules=extensions,
    install_requires=['numpy'],
    python_requires='>=3.5',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords='science astrophysics accretion',
)

