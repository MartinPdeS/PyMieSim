#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU

# Configure the medium, shell and core materials for the sphere
core_properties = [Material.silver, Material.fused_silica, 1.4 * RIU]
shell_properties = [Material.silver, Material.fused_silica, 1.4 * RIU]
medium_properties = [Material.water, 1.1 * RIU]

# Define measures to test
measures = CoreShell.available_measure_list

@pytest.mark.parametrize('medium_property', medium_properties, ids=[f'Medium:{m}' for m in medium_properties])
@pytest.mark.parametrize('core_property', core_properties, ids=[f'Property:{m}' for m in core_properties])
@pytest.mark.parametrize('shell_property', shell_properties, ids=[f'Property:{m}' for m in shell_properties])
@pytest.mark.parametrize('measure', measures)
def test_measure(measure, core_property, shell_property, medium_property):
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
        shell_property=shell_property,
        core_property=core_property,
        medium_property=medium_property
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
