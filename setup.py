#!/usr/bin/env python
# -*- coding: utf-8 -*-

# docker run -ti -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /bin/bash
# yum install wget
# wget -c 'http://sourceforge.net/projects/boost/files/boost/1.75.0/boost_1_75_0.tar.bz2'
# tar xf boost_1_75_0.tar.bz2
# cp -r boost_1_75_0/boost /usr/include
# mkdir GitProject
# cd GitProject
# git clone https://github.com/MartinPdeS/PyMieSim.git
# cd PyMieSim
# sudo python3 -m build
# mkdir /output
# /opt/python/cp37-37m/pip wheel . -w output
# auditwheel repair /output/mylibrary*whl -w /output
# python3 -m twine upload --repository pypi dist/
# python3 -m twine upload --repository pypi dist/*.tar.gz


import io
import os
import sys
import numpy
from shutil import rmtree
import pathlib
from setuptools import setup, Extension
import subprocess
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.build_ext import build_ext

from setuptools import find_packages, setup, Command
from setuptools.command.build_ext import build_ext as _build_ext


# Package meta-data.
NAME            = 'PyMieSim'
DESCRIPTION     = 'A package for light scattering simulations.'
URL             = 'https://github.com/MartinPdeS/PyMieSim'
EMAIL           = 'Martin.poinsinet.de.sivry@gmail.com'
AUTHOR          = 'Martin Poinsinet de Sivry',
REQUIRES_PYTHON = '>3.6.0'
VERSION         = '0.2.8'

# What packages are required for this module to be executed?
REQUIRED = ['scipy',
            'matplotlib',
            'cython',
            'pybind11',
            'vtk',
            'pandas',
            'beartype',
            'mayavi',
            'coverage',
            'vtk',
            'numpy']


class get_pybind11_include(object):
    """Defer numpy.get_include() until after numpy is installed."""
    def __str__(self):
        import pybind11
        return pybind11.get_include()

class get_numpy_include(object):
    """Defer numpy.get_include() until after numpy is installed."""

    def __str__(self):
        import numpy
        return numpy.get_include()



EXTRAS = {}

here = os.path.abspath(os.path.dirname(__file__))


macro = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

compile_args=["-std=c++14",
              '-fopenmp',
              '-O3',
              '-I/usr/include/boost',
              '-march=native']

link_args=['-I/usr/include/boost',]

ext_modules = [
                Extension(name               = "PyMieSim._Coupling",
                          sources            = ["PyMieSim/Coupling.pyx"],
                          include_dirs        = [numpy.get_include()],
                          language            = "c++",
                          define_macros       = macro,
                          extra_compile_args  = compile_args,
                          extra_link_args     = link_args),

                ]


try:
    with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


# Where the magic happens:
setup(
    name                          = NAME,
    version                       = about['__version__'],
    description                   = DESCRIPTION,
    long_description              = long_description,
    long_description_content_type = 'text/markdown',
    author                        = AUTHOR,
    author_email                  = EMAIL,
    setup_requires                = ['numpy', 'pybind11','cython'],
    python_requires               = '>=3.6',#REQUIRES_PYTHON,
    url                           = URL,
    packages                      = find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    install_requires              = REQUIRED,
    extras_require                = EXTRAS,
    dependency_links              = ["fibermodes @ git+https://github.com/cbrunet/fibermodes.git#egg=fibermodes-0.2.0"],
    include_package_data          = True,
    ext_modules                   = ext_modules,
    license                       = 'MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Programming Language :: C++',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
    ],
    # $ setup.py publish support.
    cmdclass={'upload': UploadCommand}, #'build_ext': build_ext
)
