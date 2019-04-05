#!/usr/bin/env python

import os

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


# https://stackoverflow.com/questions/42585210/extending-setuptools-extension-to-use-cmake-in-setup-py
class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])


class BuildExt(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = os.getcwd()

        os.makedirs(self.build_temp, exist_ok=True)
        extpath = os.path.abspath(self.get_ext_fullpath(ext.name))
        extpath_dir = os.path.dirname(extpath)
        extname, extsuffix = os.path.splitext(os.path.basename(extpath))
        os.makedirs(extpath_dir, exist_ok=True)
        target = ext.name.split('.')[-1]

        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}'.format(extpath_dir),
            '-DMODULE_OUTPUT_NAME={}'.format(extname),
            '-DMODULE_OUTPUT_SUFFIX={}'.format(extsuffix),
            '-DCMAKE_BUILD_TYPE={}'.format(config)
        ]

        build_args = [
            '--config', config,
            '--target', target,
            '--', '-j4',
        ]

        os.chdir(self.build_temp)
        self.spawn(['cmake', cwd] + cmake_args)
        self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(cwd)


setup(
    name='freddi',
    version='2.0a0.dev0',
    url='http://xray.sai.msu.ru/~malanchev/freddi/',
    license='MIT',
    author='Konstantin Malanchev',
    author_email='malanchev@sai.msu.ru',
    description='Compute FRED light curves of LMXBs outbursts',
    ext_modules=[CMakeExtension('freddi._freddi')],
    package_dir={'': 'python'},
    packages=['freddi'],
    cmdclass={'build_ext': BuildExt},
    install_requires=['numpy'],
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

