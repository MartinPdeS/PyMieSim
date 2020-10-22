#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    print('Please install or upgrade setuptools or pip to continue')
    sys.exit(1)


requirements = [ 'numpy',
                'matplotlib',
                'PyMieScatt',
                'progress',
                'tqdm',
                'fibermodes @ git+https://github.com/cbrunet/fibermodes.git#egg=0.2.0'
                ]

setup_dict = dict(
      name = 'MieCoupling',
      description='Coupled mode modlisation for fiber optic coupler',
      version = '0.1',
      author = 'Martin Poinsinet de Sivry',
      author_email = 'Martin.poinsinet.de.sivry@gmail.com',
      packages=['PyMieCoupling',
                'PyMieCoupling.classes',
                'PyMieCoupling.functions'],
      install_requires = requirements,
      license = 'Full private, no reproduction authorized',
      #url='https://gitlab.com/PolyMtlLFO/SuPyModes',
      long_description=open('README.md').read(),
      platforms = ['Linux', 'Max OSX']
)

setup(**setup_dict)
