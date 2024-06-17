#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import PyMieSim.experiment.measure as pms_measure
from PyOptik import UsualMaterial

# Define core materials and properties
core_options = [
    {'name': 'BK7', 'properties': {'core_material': UsualMaterial.BK7}},
    {'name': 'Silver', 'properties': {'core_material': UsualMaterial.Silver}},
    {'name': 'Aluminium', 'properties': {'core_material': UsualMaterial.Aluminium}},
    {'name': 'Index', 'properties': {'core_index': 1.4}}
]

# Define shell materials and properties
shell_options = [
    {'name': 'BK7', 'properties': {'shell_material': UsualMaterial.BK7}},
    {'name': 'Silver', 'properties': {'shell_material': UsualMaterial.Silver}},
    {'name': 'Aluminium', 'properties': {'shell_material': UsualMaterial.Aluminium}},
    {'name': 'Index', 'properties': {'shell_index': 1.4}}
]

# Define medium materials and properties
medium_options = [
    {'name': 'BK7', 'properties': {'medium_material': UsualMaterial.Water}},
    {'name': 'BK7', 'properties': {'medium_material': UsualMaterial.Water}},
    {'name': 'Index', 'properties': {'medium_index': 1.1}}
]

# Define measures to test
measures = pms_measure.__coreshell__


@pytest.mark.parametrize('medium_config', [m['properties'] for m in medium_options], ids=[m['name'] for m in medium_options])
@pytest.mark.parametrize('shell_config', [s['properties'] for s in shell_options], ids=[s['name'] for s in shell_options])
@pytest.mark.parametrize('core_config', [c['properties'] for c in core_options], ids=[c['name'] for c in core_options])
@pytest.mark.parametrize('measure', measures.values(), ids=measures.keys())
def test_coreshell_scattering_properties(measure, medium_config, core_config, shell_config):
    # Setup Gaussian source
    source = Gaussian(
        wavelength=np.linspace(400e-9, 1800e-9, 50),
        polarization=0,
        optical_power=1e-3,
        NA=0.2
    )

    # Setup core-shell scatterer
    scatterer = CoreShell(
        core_diameter=np.linspace(400e-9, 1400e-9, 10),
        shell_width=300e-9,
        source=source,
        **medium_config,
        **core_config,
        **shell_config
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

    experiment.get(measure)


if __name__ == "__main__":
    pytest.main()


# -
