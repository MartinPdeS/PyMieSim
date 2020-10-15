#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    print('Please install or upgrade setuptools or pip to continue')
    sys.exit(1)


requirements = [ 'numpy', 'matplotlib', 'PyMieScatt', 'progress']

setup_dict = dict(
      description='Coupled mode modlisation for fiber optic coupler',
      name = 'MieCoupling',
      version = '0.1',
      author = 'Martin Poinsinet de Sivry',
      author_email = 'Martin.poinsinet.de.sivry@gmail.com',
      packages=find_packages(),
      py_modules = [],
      install_requires = requirements,
      dependency_links=['https://github.com/cbrunet/fibermodes.git#egg=package-1.0'],
      license = 'Full private, no reproduction authorized',
      #url='https://gitlab.com/PolyMtlLFO/SuPyModes',
      long_description=open('README.md').read(),
      platforms = ['Linux', 'Max OSX']
)

setup(**setup_dict)
