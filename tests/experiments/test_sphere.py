#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.experiment.measure as pms_measure
from PyOptik import materials

# Configure the core materials for the sphere
core_options = [
    {'name': 'crown', 'properties': {'material': materials.crown}},
    {'name': 'fused silica', 'properties': {'material': materials.fused_silica}},
    {'name': 'Index', 'properties': {'index': 1.4}}
]

# Define medium options
medium_options = [
    {'name': 'water', 'properties': {'medium_material': materials.water}},
    {'name': 'Index', 'properties': {'medium_index': 1.1}}
]

# List of measures to be tested
measures = pms_measure.__sphere__


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures.values(), ids=measures.keys())
def test_sphere_scattering_properties(measure, core_config, medium_config):
    # Set up the Gaussian source
    source = Gaussian(
        wavelength=np.linspace(600e-9, 1000e-9, 50),
        polarization=0,
        optical_power=1e-3,
        NA=0.2
    )

    # Configure the spherical scatterer
    scatterer = Sphere(
        diameter=np.linspace(400e-9, 1400e-9, 10),
        source=source,
        **medium_config,
        **core_config
    )

    # Configure the detector
    detector = CoherentMode(
        mode_number='LP01',
        rotation=0,
        NA=0.2,
        polarization_filter=None,
        gamma_offset=0,
        phi_offset=0,
        sampling=100
    )

    # Set up and run the experiment
    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    experiment.get(measure)


if __name__ == "__main__":
    pytest.main([__file__])

# -
