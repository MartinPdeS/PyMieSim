#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.measure as pms_measure
from PyMieSim.materials import Silver, BK7, Aluminium

# Configure the core materials for the sphere
core_options = [
    {'name': 'BK7', 'properties': {'material': BK7}},
    {'name': 'Silver', 'properties': {'material': Silver}},
    {'name': 'Aluminium', 'properties': {'material': Aluminium}},
    {'name': 'Index', 'properties': {'index': 1.4}}
]

# Define medium options
medium_options = [
    {'name': 'BK7', 'properties': {'medium_material': BK7}},
    {'name': 'Index', 'properties': {'medium_index': 1.1}}
]

# List of measures to be tested
measures = [
    pms_measure.Qsca, pms_measure.Qabs, pms_measure.Qback,
    pms_measure.g, pms_measure.a1, pms_measure.b1, pms_measure.coupling
]


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures, ids=[m.short_label for m in measures])
def test_sphere_scattering_properties(measure, core_config, medium_config):
    # Set up the Gaussian source
    source_set = Gaussian(
        wavelength=np.linspace(400e-9, 1800e-9, 50),
        polarization_value=0,
        polarization_type='linear',
        optical_power=1e-3,
        NA=0.2
    )

    # Configure the spherical scatterer
    scatterer_set = Sphere(
        diameter=np.linspace(400e-9, 1400e-9, 10),
        source_set=source_set,
        **medium_config,
        **core_config
    )

    # Configure the detector
    detector_set = Photodiode(
        NA=0.2,
        polarization_filter=None,
        gamma_offset=0,
        phi_offset=0,
        sampling=100
    )

    # Set up and run the experiment
    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=detector_set
    )

    # Perform the measurement
    experiment.get(measure)

# -
