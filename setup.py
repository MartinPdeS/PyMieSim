#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
from setuptools import setup, find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

requirements = ['numpy',
                #'cupy',
                'matplotlib',
                'cython',
                'pandas',
                'tqdm',
                'fibermodes @ git+https://github.com/cbrunet/fibermodes.git#egg=0.2.0',
                ]


ext_modules = [ Extension("PyMieCoupling.cython.S1S2",
                         ["PyMieCoupling/cython/S1S2.pyx"],
                         include_dirs = ['.'],
                         language="c++",
                         extra_compile_args=["-std=c++11", '-fopenmp', '-lboost_filesystem','-lboost_system'],
                         extra_link_args=["-std=c++11",'-fopenmp','-lboost_filesystem','-lboost_system']),

                Extension("PyMieCoupling.cython.Fields",
                         ["PyMieCoupling/cython/Fields.pyx"],
                         include_dirs = ['.'],
                         language="c++",
                         extra_compile_args=["-std=c++11", '-fopenmp','-lboost_filesystem','-lboost_system'],
                         extra_link_args=["-std=c++11",'-fopenmp','-lboost_filesystem','-lboost_system']),
                         ]

setup_dict = dict(
      name               = 'PyMieCoupling',
      description        = 'Coupled mode modlisation for fiber optic coupler',
      version            = '0.1',
      author             = 'Martin Poinsinet de Sivry',
      author_email       = 'Martin.poinsinet.de.sivry@gmail.com',
      cmdclass           = {'build_ext': build_ext},
      packages           = ['PyMieCoupling',
                          'PyMieCoupling.classes',
                          'PyMieCoupling.functions'],
      install_requires   = requirements,
      license            = 'Full private, no reproduction authorized',
      url                = 'https://gitlab.com/PolyMtlLFO/SuPyModes',
      long_description   = open('README.md').read(),
      platforms          = ['Linux', 'Max OSX'],
      include_dirs       = [numpy.get_include()],
      ext_modules        = ext_modules,
      zip_safe           = False)


setup(**setup_dict)
