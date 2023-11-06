#!/usr/bin/env python
# -*- coding: utf-8 -*-

from DataVisual import Xparameter


g = Xparameter(name='g', format="<20s", unit="1", long_label='Anisotropy coefficient', short_label='g')
coupling = Xparameter(name='coupling', format="<20s", unit="Watt", long_label='Coupling', short_label='Coupling')

Qsca = Xparameter(name='Qsca', format="<20s", unit="1", long_label='Scattering efficiency', short_label='Qsca')
Qext = Xparameter(name='Qext', format="<20s", unit="1", long_label='Extinction efficiency', short_label='Qext')
Qabs = Xparameter(name='Qabs', format="<20s", unit="1", long_label='Absorption efficiency', short_label='Qabs')
Qratio = Xparameter(name='Qratio', format="<20s", unit="1", long_label='Ratio efficiency', short_label='Qratio')
Qback = Xparameter(name='Qback', format="<20s", unit="1", long_label='Backward efficiency', short_label='Qback')
Qforw = Xparameter(name='Qforw', format="<20s", unit="1", long_label='Forward efficiency', short_label='Qforw')
Qpr = Xparameter(name='Qpr', format="<20s", unit="1", long_label='Radiation pressure efficiency', short_label='Qpr')

Csca = Xparameter(name='Csca', format="<20s", unit="m^2", long_label='Scattering cross-section', short_label='Csca')
Cext = Xparameter(name='Cext', format="<20s", unit="m^2", long_label='Extinction cross-section', short_label='Cext')
Cabs = Xparameter(name='Cabs', format="<20s", unit="m^2", long_label='Absorption cross-section', short_label='Cabs')
Cratio = Xparameter(name='Cratio', format="<20s", unit="m^2", long_label='Ratio cross-section', short_label='Cratio')
Cback = Xparameter(name='Cback', format="<20s", unit="m^2", long_label='Backward cross-section', short_label='Cback')
Cforw = Xparameter(name='Cforw', format="<20s", unit="m^2", long_label='Forward cross-section', short_label='Cforw')
Cpr = Xparameter(name='Cpr', format="<20s", unit="m^2", long_label='Radiation pressure cross-section', short_label='Cpr')

a1 = Xparameter(name='a1', format="<20s", unit=None, long_label='Electric dipole coefficient')
a2 = Xparameter(name='a2', format="<20s", unit=None, long_label='Electric quadrupole coefficient')
a3 = Xparameter(name='a3', format="<20s", unit=None, long_label='Electric octopole coeffcient')
b1 = Xparameter(name='b1', format="<20s", unit=None, long_label='Magnertic dipole coefficient')
b2 = Xparameter(name='b2', format="<20s", unit=None, long_label='Magnetic quadrupole coefficient')
b3 = Xparameter(name='b3', format="<20s", unit=None, long_label='Magnetic octopole coefficient')

a11 = Xparameter(name='a11', format="<20s", unit=None, long_label='Electric dipole coefficient')
a21 = Xparameter(name='a21', format="<20s", unit=None, long_label='Electric quadrupole coefficient')
a12 = Xparameter(name='a12', format="<20s", unit=None, long_label='Electric octopole coeffcient')
a22 = Xparameter(name='a22', format="<20s", unit=None, long_label='Electric dipole coefficient')
a13 = Xparameter(name='a13', format="<20s", unit=None, long_label='Electric quadrupole coefficient')
a23 = Xparameter(name='a23', format="<20s", unit=None, long_label='Electric octopole coeffcient')

b11 = Xparameter(name='b11', format="<20s", unit=None, long_label='Electric dipole coefficient')
b21 = Xparameter(name='b21', format="<20s", unit=None, long_label='Electric quadrupole coefficient')
b12 = Xparameter(name='b12', format="<20s", unit=None, long_label='Electric octopole coeffcient')
b22 = Xparameter(name='b22', format="<20s", unit=None, long_label='Electric dipole coefficient')
b13 = Xparameter(name='b13', format="<20s", unit=None, long_label='Electric quadrupole coefficient')
b23 = Xparameter(name='b23', format="<20s", unit=None, long_label='Electric octopole coeffcient')

# -
