#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.units import nanometer, degree, watt, AU, RIU, volt, meter

# Configure the medium and core materials for the sphere
properties = [Material.silver, Material.fused_silica, 1.4 * RIU]
medium_properties = [Material.water, 1.1 * RIU]

# Measures to be tested
measures = Cylinder.available_measure_list

gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 50) * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 50) * nanometer,
    polarization=0 * degree,
    amplitude=1 * volt / meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize('medium_property', medium_properties, ids=[f'Medium:{m}' for m in medium_properties])
@pytest.mark.parametrize('source', sources, ids=[f'Source:{m.__class__.__name__}' for m in sources])
@pytest.mark.parametrize('property', properties, ids=[f'Property:{m}' for m in properties])
@pytest.mark.parametrize('measure', measures)
def test_measure(measure, source, medium_property, property):
    # Configure the spherical scatterer
    scatterer = Cylinder(
        diameter=np.linspace(400, 1400, 10) * nanometer,
        source=source,
        medium_property=medium_property,
        property=property
    )

    # Configure the detector
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
    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
