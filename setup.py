#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import sys
from setuptools import setup, find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy
from pybind11.setup_helpers import Pybind11Extension


macro = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
compile_args=["-std=c++11",
              '-fopenmp',
              '-lboost_filesystem',
              '-lboost_system',
              '-O3',
              '-march=native']

link_args=["-std=c++11",
           '-fopenmp',
           '-lboost_filesystem',
           '-lboost_system',
           '-O3',
           '-march=native']

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


ext_modules = [ Extension(name               = "PyMieSim.cython.S1S2",
                          sources            = ["PyMieSim/cython/S1S2.pyx"],
                          include_dirs       = ['.'],
                          language           = "c++",
                          define_macros      = macro,
                          extra_compile_args = compile_args,
                          extra_link_args    = link_args
                         ),

                Extension(name               = "PyMieSim.cpp.Interface",
                          sources            = ["PyMieSim/cpp/Interface.pxd"],
                          include_dirs       = [numpy.get_include()],
                          language           = "c++",
                          define_macros      = macro,
                          extra_compile_args = compile_args,
                          extra_link_args    = link_args
                         ),


                Extension(name               = "PyMieSim.cpp.Sphere",
                          sources            = ["PyMieSim/cpp/Sphere.pyx"],
                          include_dirs       = [numpy.get_include()],
                          language           = "c++",
                          define_macros      = macro,
                          extra_compile_args = compile_args,
                          extra_link_args    = link_args
                         ),




                Extension(name               = "PyMieSim.cpp.Coupling",
                          sources            = ["PyMieSim/cpp/Coupling.pyx"],
                          include_dirs        = [numpy.get_include()],
                          language            = "c++",
                          define_macros       = macro,
                          extra_compile_args  = compile_args,
                          extra_link_args     = link_args
                         ),


                Pybind11Extension(name            = "GLMT",
                                  sources         = ["PyMieSim/cpp/GLMT.cpp"],
                                  language        = 'c++11',
                                  extra_link_args = link_args
                ),


                         ]

setup_dict = dict(
      name               = 'PyMieSim',
      description        = 'A package for light scattering simulations',
      version            = '0.1',
      #python_requires    = ">=3.0*",
      install_requires   = ['numpy',
                            'scipy',
                            'matplotlib',
                            'pandas',
                            'cython',
                            'mayavi',
                            'ai.cs',
                            'geos',
                            'cartopy',
                            "fibermodes @ git+https://github.com/cbrunet/fibermodes.git#egg=fibermodes-0.2.0",
                            ],
      #dependency_links   = ['https://github.com/cbrunet/fibermodes.git#egg=fibermodes-0.2.0'],
      author             = 'Martin Poinsinet de Sivry',
      author_email       = 'Martin.poinsinet.de.sivry@gmail.com',
      cmdclass           = {'build_ext': build_ext},
      packages           = ['PyMieSim',
                           'PyMieSim.functions'],
      license            = 'Full private, no reproduction authorized',
      url                = 'https://gitlab.com/Martth/miecoupling',
      long_description   = open('README.md').read(),
      platforms          = ['Linux', 'Max OSX'],
      include_dirs       = [numpy.get_include()],
      ext_modules        = ext_modules,
      zip_safe           = False)


setup(**setup_dict)
