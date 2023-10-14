#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Distribution
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


import os

try:
    os.add_dll_directory('C:/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/lib')
except Exception:
    pass


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


class bdist_wheelNotPure(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        self.root_is_pure = False


setup(cmdclass={'bdist_wheelPure': _bdist_wheel, 'bdist_wheel': bdist_wheelNotPure},
      distclass=BinaryDistribution,
      )
