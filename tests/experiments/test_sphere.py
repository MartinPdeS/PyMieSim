#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Configure the core materials for the sphere
core_options = [
    {'name': 'crown', 'properties': {'material': Material.silver}},
    {'name': 'fused silica', 'properties': {'material': Material.fused_silica}},
    {'name': 'Index', 'properties': {'index': 1.4 * RIU}}
]

# Define medium options
medium_options = [
    {'name': 'water', 'properties': {'medium_material': Material.water}},
    {'name': 'Index', 'properties': {'medium_index': 1.1 * RIU}}
]

# List of measures to be tested
measures = [
    'Qsca', 'Qext', 'Qabs', 'Qback', 'Qforward', 'Qratio',
    'g', 'Qpr', 'Csca', 'Cext', 'Cabs', 'Cratio',
    'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'coupling'
]

@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures)
def test_sphere_scattering_properties(measure, core_config, medium_config):
    # Set up the Gaussian source
    source = Gaussian(
        wavelength=np.linspace(600, 1000, 50) * nanometer,
        polarization=0 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    # Configure the spherical scatterer
    scatterer = Sphere(
        diameter=np.linspace(400, 1400, 10) * nanometer,
        source=source,
        **medium_config,
        **core_config
    )

    # Configure the detector
    detector = CoherentMode(
        mode_number='LP01',
        rotation=0 * degree,
        NA=0.2 * AU,
        polarization_filter=None,
        gamma_offset=0 * degree,
        phi_offset=0 * degree,
        sampling=100 * AU
    )

    # Set up and run the experiment
    experiment = Setup(scatterer=scatterer, source=source, detector=detector)

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main([__file__])

# -
