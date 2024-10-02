#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.experiment.measure as pms_measure
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Define core materials and properties
core_options = [
    {'name': 'BK7', 'properties': {'core_material': Material.BK7}},
    {'name': 'Index', 'properties': {'core_index': 1.4 * RIU}}
]

# Define shell materials and properties
shell_options = [
    {'name': 'BK7', 'properties': {'shell_material': Material.BK7}},
    {'name': 'Index', 'properties': {'shell_index': 1.4 * RIU}}
]

# Define medium materials and properties
medium_options = [
    {'name': 'water', 'properties': {'medium_material': Material.water}},
    {'name': 'Index', 'properties': {'medium_index': 1.1 * RIU}}
]

# Define measures to test
measures = [
    'Qsca', 'Qext', 'Qabs', 'Qback', 'Qforward', 'Qratio',
    'g', 'Qpr', 'Csca', 'Cext', 'Cabs', 'Cratio',
    'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'coupling'
]


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('shell_config', [s['properties'] for s in shell_options], ids=[s['name'] for s in shell_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures)
def test_coreshell_scattering_properties(measure, medium_config, core_config, shell_config):
    # Setup Gaussian source
    source = Gaussian(
        wavelength=np.linspace(800, 1000, 50) * nanometer,
        polarization=0 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    # Setup core-shell scatterer
    scatterer = CoreShell(
        core_diameter=np.linspace(800, 1000, 10) * nanometer,
        shell_width=300 * nanometer,
        source=source,
        **medium_config,
        **core_config,
        **shell_config
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

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main([__file__])


# -
