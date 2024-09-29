#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.experiment.measure as pms_measure
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Material configurations for the cylinder core
core_options = [
    {'name': 'BK7', 'properties': {'material': Material.BK7}},
    {'name': 'Index', 'properties': {'index': 1.4 * RIU}}
]

# Medium configurations
medium_options = [
    {'name': 'water', 'properties': {'medium_material': Material.water}},
    {'name': 'Index', 'properties': {'medium_index': 1.1 * RIU}}
]

# Measures to be tested
measures = pms_measure.__cylinder__


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures.values(), ids=measures.keys())
def test_cylinder_scattering_properties(measure, medium_config, core_config):
    # Setup Gaussian source
    source = Gaussian(
        wavelength=np.linspace(600, 1000, 50) * nanometer,
        polarization=0 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    # Setup cylindrical scatterer
    scatterer = Cylinder(
        diameter=np.linspace(400, 1400, 10) * nanometer,
        source=source,
        **medium_config,
        **core_config
    )

    # Setup detector
    detector = Photodiode(
        NA=0.2 * AU,
        polarization_filter=None,
        gamma_offset=0 * degree,
        phi_offset=0 * degree,
        sampling=100 * AU
    )

    # Configure and run the experiment
    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    # Execute measurement
    experiment.get(measure)


if __name__ == "__main__":
    pytest.main([__file__])

# -
