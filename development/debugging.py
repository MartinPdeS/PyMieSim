#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import matplotlib.pyplot as plt
from PyMieSim import single
from PyMieSim import experiment, measure

source = single.source.Gaussian(
    wavelength=1e-6,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.3
)

scatterer = single.scatterer.Sphere(
    diameter=1e-6,
    index=1.5,
    source=source
)

kwargs = dict(
    mode_number='LP21',
    NA=0.1,
    gamma_offset=0,
    phi_offset=0,
    polarization_filter=None,
)

array = numpy.zeros(100)
angle_degree = numpy.linspace(0, 180, 100)
angle_radian = numpy.linspace(0, numpy.pi, 100)
for idx, rotation in enumerate(angle_degree):
    detector = single.detector.CoherentMode(**kwargs, rotation=rotation)
    array[idx] = detector.coupling(scatterer)

figure, ax = plt.subplots(1, 1)

ax.plot(array)
plt.show()