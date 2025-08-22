#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment import Setup
from PyOptik import Material

# Configure the medium and core materials for the sphere
properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]

# List of measures to be tested
measures = Sphere.available_measure_list

gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 150) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 150) * ureg.nanometer,
    polarization=0 * ureg.degree,
    amplitude=1 * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize('medium_property', medium_properties, ids=[f'Medium:{m}' for m in medium_properties])
@pytest.mark.parametrize('source', sources, ids=[f'Source:{m.__class__.__name__}' for m in sources])
@pytest.mark.parametrize('property', properties, ids=[f'Property:{m}' for m in properties])
@pytest.mark.parametrize('measure', measures)
def test_get_measure(source, measure, property, medium_property):
    # Configure the spherical scatterer
    scatterer = Sphere(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        source=source,
        property=property,
        medium_property=medium_property
    )

    # Configure the detector
    detector = CoherentMode(
        mode_number='LP01',
        rotation=0 * ureg.degree,
        NA=[0.1] * ureg.AU,
        polarization_filter=None,
        gamma_offset=[0, 1] * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100 * ureg.AU
    )

    # Set up and run the experiment
    experiment = Setup(scatterer=scatterer, source=source, detector=detector)

    experiment.get(measure, drop_unique_level=True, scale_unit=True)

    experiment.get(measure, drop_unique_level=False, scale_unit=True, as_numpy=True)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
