#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyOptik import Material
from TypedUnit import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment import Setup

# Configure the medium, shell and core materials for the sphere
core_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
shell_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]

# Define measures to test
measures = CoreShell.available_measure_list

gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    amplitude=1 * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize('medium_property', medium_properties, ids=[f'Medium:{m}' for m in medium_properties])
@pytest.mark.parametrize('core_property', core_properties, ids=[f'Property:{m}' for m in core_properties])
@pytest.mark.parametrize('source', sources, ids=[f'Source:{m.__class__.__name__}' for m in sources])
@pytest.mark.parametrize('shell_property', shell_properties, ids=[f'Property:{m}' for m in shell_properties])
@pytest.mark.parametrize('measure', measures)
def test_measure(measure, source, core_property, shell_property, medium_property):
    # Setup core-shell scatterer
    scatterer = CoreShell(
        core_diameter=np.linspace(800, 1000, 10) * ureg.nanometer,
        shell_thickness=300 * ureg.nanometer,
        source=source,
        shell_property=shell_property,
        core_property=core_property,
        medium_property=medium_property
    )

    # Setup detector
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        polarization_filter=None,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100 * ureg.AU
    )

    # Configure and run the experiment
    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
