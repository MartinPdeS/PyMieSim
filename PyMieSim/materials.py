#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik import DataMeasurement, Sellmeier

BK7 = Sellmeier('BK7')
FusedSilica = Sellmeier('silica')
SodaLimeGlass = DataMeasurement('sodalimeglass')
Silver = DataMeasurement('silver')
Gold = DataMeasurement('gold')
Aluminium = DataMeasurement('aluminium')
SI = Sellmeier('silica')
SIO2 = DataMeasurement('sio2')
TIO2 = DataMeasurement('tio2')
