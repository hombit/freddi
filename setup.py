#!/usr/bin/env python

from skbuild import setup


with open('Readme.md') as fh:
    readme = fh.read()


setup(
    name='freddi',
    version='2.0.0b2',
    url='http://xray.sai.msu.ru/~malanchev/freddi/',
    license='GPLv3',
    author='Konstantin Malanchev',
    author_email='malanchev@sai.msu.ru',
    description='Compute FRED light curves of LMXBs outbursts',
    long_description=readme,
    long_description_content_type='text/markdown',
    package_dir={'': 'python'},
    packages=['freddi'],
    install_requires=['numpy'],
    setup_requires=['numpy'],
    python_requires='>=3.7',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords='science astrophysics accretion',
)

