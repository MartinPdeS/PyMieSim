#!/usr/bin/env python
# -*- coding: utf-8 -*-

from DataVisual import units

# Efficiencies
Qsca = units.Efficiency(long_label='Scattering Efficiencies', short_label='Qsca')
Qext = units.Efficiency(long_label='Extinction Efficiencies', short_label='Qext')
Qabs = units.Efficiency(long_label='Absorption Efficiencies', short_label='Qabs')
Qratio = units.Efficiency(long_label='Ratio Efficiencies', short_label='Qratio')
Qforw = units.Efficiency(long_label='Forward Efficiencies', short_label='Qforw')
Qback = units.Efficiency(long_label='Backward Efficiencies', short_label='Qback')
Qpr = units.Efficiency(long_label='Radiation pres Efficiencies', short_label='Qpr')

# Cross section
Csca = units.Area(long_label='Scattering Cross-Section', short_label='Csca')
Cext = units.Area(long_label='Extinction Cross-Section', short_label='Cext')
Cabs = units.Area(long_label='Absorption Cross-Section', short_label='Cabs')
Cratio = units.Area(long_label='Ratio Cross-Section', short_label='Cratio')
Cforw = units.Area(long_label='Forward Cross-Section', short_label='Cforw')
Cback = units.Area(long_label='Backward Cross-Section', short_label='Cback')
Cpr = units.Area(long_label='Radiation pres Cross-Section', short_label='Cpr')

# Electric & magnetic multipole
a1 = units.Custom(short_label='a1', long_label='Electric dipole coefficient')
a2 = units.Custom(short_label='a2', long_label='Electric quadrupole coefficient')
a3 = units.Custom(short_label='a3', long_label='Electric octopole coeffcient')
b1 = units.Custom(short_label='b1', long_label='Magnetic dipole coefficient')
b2 = units.Custom(short_label='b2', long_label='Magnetic quadrupole coefficient')
b3 = units.Custom(short_label='b3', long_label='Magnetic octopole coefficient')

# Cylinder Electric & magnetic multipole
a11 = units.Custom(short_label='a11', long_label='Electric dipole coefficient')
a21 = units.Custom(short_label='a21', long_label='Electric quadrupole coefficient')
a12 = units.Custom(short_label='a12', long_label='Electric octopole coeffcient')
a22 = units.Custom(short_label='a22', long_label='Electric dipole coefficient')
a13 = units.Custom(short_label='a13', long_label='Electric quadrupole coefficient')
a23 = units.Custom(short_label='a23', long_label='Electric octopole coeffcient')

b11 = units.Custom(short_label='b11', long_label='Magnetic dipole coefficient')
b21 = units.Custom(short_label='b21', long_label='Magnetic quadrupole coefficient')
b12 = units.Custom(short_label='b12', long_label='Magnetic octopole coeffcient')
b22 = units.Custom(short_label='b22', long_label='Magnetic dipole coefficient')
b13 = units.Custom(short_label='b13', long_label='Magnetic quadrupole coefficient')
b23 = units.Custom(short_label='b23', long_label='Magnetic octopole coeffcient')


# Extra
g = units.Custom(short_label='g', long_label='Anisotropy coefficient')
coupling = units.Power(short_label='coupling', long_label='Coupling')

# -
