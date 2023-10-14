#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import sqrt, pi

c = 299792458.0  #: Speed of light in vacuum (m/s).
h = 6.62606957e-34  #: Plank constant (m²kg/s).
mu0 = 1.2566370614359173e-06  #: Vacuum permeability (H/m).
epsilon0 = 8.854187817620389e-12  #: Vacuum permittivity (F/m).
eV = 1.602176565e-19  #: Electron charge (C).

tpi = 2 * pi  #: Two times pi
eta0 = sqrt(mu0 / epsilon0)  #: Impedance of free-space.
Y0 = sqrt(epsilon0 / mu0)  #: Admitance of free-space.
