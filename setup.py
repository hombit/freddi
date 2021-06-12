#!/usr/bin/env python

import os

from skbuild import setup
from setuptools.command.build_ext import build_ext


with open('Readme.md') as fh:
    readme = fh.read()


setup(
    name='freddi',
    version='2.0b0',
    url='http://xray.sai.msu.ru/~malanchev/freddi/',
    license='MIT',
    author='Konstantin Malanchev',
    author_email='malanchev@sai.msu.ru',
    description='Compute FRED light curves of LMXBs outbursts',
    long_description=readme,
    long_description_content_type='text/markdown',
    package_dir={'': 'python'},
    packages=['freddi'],
    install_requires=['numpy'],
    setup_requires=['numpy'],
    tests_require=['parameterized'],
    python_requires='>=3.5',
    test_suite='test',
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

