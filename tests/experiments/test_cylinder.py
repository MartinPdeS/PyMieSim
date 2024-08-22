#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.experiment.measure as pms_measure
from PyOptik import UsualMaterial

# Material configurations for the cylinder core
core_options = [
    {'name': 'Silver', 'properties': {'material': UsualMaterial.Silver}},
    {'name': 'Aluminium', 'properties': {'material': UsualMaterial.Aluminium}},
    {'name': 'Index', 'properties': {'index': 1.4}}
]

# Medium configurations
medium_options = [
    {'name': 'BK7', 'properties': {'medium_material': UsualMaterial.BK7}},
    {'name': 'Index', 'properties': {'medium_index': 1.1}}
]

# Measures to be tested
measures = pms_measure.__cylinder__


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures.values(), ids=measures.keys())
def test_cylinder_scattering_properties(measure, medium_config, core_config):
    # Setup Gaussian source
    source = Gaussian(
        wavelength=np.linspace(400e-9, 1800e-9, 50),
        polarization=0,
        optical_power=1e-3,
        NA=0.2
    )

    # Setup cylindrical scatterer
    scatterer = Cylinder(
        diameter=np.linspace(400e-9, 1400e-9, 10),
        source=source,
        **medium_config,
        **core_config
    )

    # Setup detector
    detector = Photodiode(
        NA=0.2,
        polarization_filter=None,
        gamma_offset=0,
        phi_offset=0,
        sampling=100
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
    pytest.main()

# -
