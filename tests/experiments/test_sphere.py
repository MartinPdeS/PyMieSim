#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy
import PyMieSim

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

from PyMieSim.materials import Silver, BK7, Aluminium


core_type = [
    {'name': 'BK7', 'kwarg': {'material': BK7}},
    {'name': 'Silver', 'kwarg': {'material': Silver}},
    {'name': 'Aluminium', 'kwarg': {'material': Aluminium}},
    {'name': 'Index', 'kwarg': {'index': 1.4}}
]

measures = [
    PyMieSim.measure.Qsca,
    PyMieSim.measure.Qabs,
    PyMieSim.measure.Qback,
    PyMieSim.measure.g,
    PyMieSim.measure.a1,
    PyMieSim.measure.b1,
    PyMieSim.measure.coupling
]


@pytest.mark.parametrize('core_type', [p['kwarg'] for p in core_type], ids=[p['name'] for p in core_type])
@pytest.mark.parametrize('measure', measures, ids=[p.name for p in measures])
def test_sphere_experiment(measure, core_type):
    source_set = Gaussian(
        wavelength=numpy.linspace(400e-9, 1800e-9, 50),
        polarization_value=0,
        polarization_type='linear',
        optical_power=1e-3,
        NA=0.2
    )

    scatterer_set = Sphere(
        n_medium=1,
        diameter=numpy.linspace(400e-9, 1400e-9, 10),
        **core_type,
        source_set=source_set
    )

    detector_set = Photodiode(
        NA=0.2,
        polarization_filter=None,
        gamma_offset=0,
        phi_offset=0,
        sampling=100
    )

    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=detector_set
    )

    experiment.get(measure)

# -
